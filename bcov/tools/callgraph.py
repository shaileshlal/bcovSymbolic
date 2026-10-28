"""Regenerate the Graphviz call graphs and the generated part of ARCHITECTURE.md.

Two independent sources of edges are merged:

* static: an ``ast`` resolver over the live modules (calls to module functions,
  imported names, ``module.attr``, ``self.method``, ``Class.method``,
  ``Class(...)`` -> ``__init__``/``__post_init__``, and dispatch tables such as
  ``quintic_frames._CONVERSIONS``). It never guesses from bare method names.
* runtime: the caller->callee pairs observed by a profiler while the test suite
  and a few CLI commands run in-process. This adds operator overloads and calls
  on objects that static resolution cannot type (``a + b``, ``m.BCOV()``).

Edges seen by both are solid, runtime-only dashed blue, static-only (real code
paths the traced runs do not reach) dotted red. Graphviz ``dot`` must be on PATH.
Run from the repository root (needs pytest; takes well under a minute):

    python bcov/tools/callgraph.py            # writes docs/callgraph/, updates ARCHITECTURE.md
    python bcov/tools/callgraph.py --no-trace # static edges only, no tests run
"""
from __future__ import annotations

import argparse
import ast
import contextlib
import io
from pathlib import Path
import subprocess
import sys
from collections.abc import Callable, Iterator

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"docs/callgraph"
DOC = ROOT/"ARCHITECTURE.md"
# Maintained public modules; historical scripts appear only in the import graph.
LIVE = ("bcov.bcov_yy", "bcov.quintic_frames", "bcov.datasets", "bcov.patterns",
        "bcov.verify_physical_qseries", "bcov.cli",
        "quintic_bcov_lab.arithmetic", "quintic_bcov_lab.sections", "quintic_bcov_lab.backend",
        "quintic_bcov_lab.datasets", "quintic_bcov_lab.runner", "quintic_bcov_lab.__main__")
SCRIPTS = ("bcov.audit_coefficients", "bcov.audit_denominators", "bcov.bcov_reconcile_frames",
           "bcov.bcov_sparse_converter", "bcov.quintic_invariant_extract", "bcov.test_frames",
           "bcov.test_gu_gauge", "bcov.verify_yy_qseries_legacy")
ENTRY_POINTS = ("Quintic_BCOV_Tutorial.py", "notebooks/arithmetic_patterns.py",
                "notebooks/physical_verification.py", "notebooks/research_overview.py",
                "run_all_checks.py", "tests/test_package.py", "tests/test_research_exact.py",
                "tests/test_research_integration.py")
CLI_RUNS = (("bcov.cli", ["recursion", "--genus", "4"]), ("bcov.cli", ["patterns", "--genus", "6"]),
            ("bcov.cli", ["verify", "--genus", "4", "--degree", "3"]),
            ("quintic_bcov_lab.__main__", ["backend-check"]), ("quintic_bcov_lab.__main__", ["list"]),
            ("quintic_bcov_lab.__main__", ["audit"]))
Edge = tuple[str, str]
EDGE_STYLES = {"both": 'style=solid', "runtime": 'style=dashed, color="#1f77b4"',
               "static": 'style=dotted, color="#d62728", penwidth=1.6'}
CLUSTER_COLOURS = ("#eef3fb", "#fdf1e6", "#eaf6ec", "#f6ecf6", "#fbf8e3", "#e8f4f6", "#f3f0ea")
POLY_HELPERS = {f"bcov.bcov_yy.{n}" for n in ("rational", "variable_index", "_weight", "_product")}
# Arithmetic/construction used everywhere; hidden as call targets outside the Poly diagram.
PRIMITIVE = POLY_HELPERS | {f"bcov.bcov_yy.Poly.{m}" for m in (
    "__init__ __add__ __sub__ __rsub__ __neg__ __mul__ __truediv__ __pow__ __eq__ "
    "__hash__ __bool__ __len__ coerce constant scale variable terms").split()}


def module_path(module: str) -> Path:
    return ROOT/(module.replace(".", "/")+".py")


# ---------------------------------------------------------------- static edges

class Module:
    """Definitions, imports and dispatch tables of one source file."""

    def __init__(self, name: str) -> None:
        self.name = name
        self.tree = ast.parse(module_path(name).read_text(encoding="utf-8"))
        self.functions: set[str] = set()                # qualified names
        self.classes: dict[str, set[str]] = {}          # class -> method names
        self.imports: dict[str, str] = {}               # local name -> module or module.symbol
        self.tables: dict[str, set[str]] = {}           # dict/tuple of functions
        package = name.rsplit(".", 1)[0]
        for node in self.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                self.functions.add(f"{name}.{node.name}")
            elif isinstance(node, ast.ClassDef):
                self.classes[node.name] = {m.name for m in node.body
                                           if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef))}
                self.functions |= {f"{name}.{node.name}.{m}" for m in self.classes[node.name]}
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ImportFrom):
                # "from .x import a", "from . import x" and the script-mode "from x import a"
                base = (package+("."+node.module if node.module else "") if node.level
                        else node.module or "")
                if not node.level and f"{package}.{base}" in LIVE+SCRIPTS:
                    base = f"{package}.{base}"
                for alias in node.names:
                    self.imports[alias.asname or alias.name] = f"{base}.{alias.name}"
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    self.imports[alias.asname or alias.name.split(".")[0]] = alias.name
        for node in self.tree.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)) and isinstance(node.value, (ast.Dict, ast.Tuple, ast.List)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                names = {f"{name}.{n.id}" for n in ast.walk(node.value)
                         if isinstance(n, ast.Name) and f"{name}.{n.id}" in self.functions}
                for t in targets:
                    if isinstance(t, ast.Name) and names:
                        self.tables[t.id] = names


def _constructor(modules: dict[str, Module], qualified: str) -> str | None:
    """__init__ (or dataclass __post_init__) of a class, or the function itself."""
    module, _, last = qualified.rpartition(".")
    m = modules.get(module)
    if m is None:
        return None
    if last in m.classes:
        for hook in ("__init__", "__post_init__"):
            if hook in m.classes[last]:
                return f"{qualified}.{hook}"
        return None
    return qualified if qualified in m.functions else None


def _resolve(modules: dict[str, Module], m: Module, cls: str | None, func: ast.expr) -> set[str]:
    """Targets of one call expression; empty when the receiver cannot be typed."""
    def symbol(name: str) -> str | None:
        if name in m.classes or f"{m.name}.{name}" in m.functions:
            return f"{m.name}.{name}"
        target = m.imports.get(name)
        if target is None:
            return None
        module, _, last = target.rpartition(".")
        if target in modules:
            return target
        if module in modules:
            imported = modules[module].imports.get(last)   # re-exported through a package
            return f"{module}.{last}" if imported is None else imported
        return None

    if isinstance(func, ast.Name):
        q = symbol(func.id)
        return {c} if q and (c := _constructor(modules, q)) else set()
    if isinstance(func, ast.Subscript) and isinstance(func.value, ast.Name):
        return set(m.tables.get(func.value.id, ()))
    if not (isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name)):
        return set()
    receiver, attr = func.value.id, func.attr
    if receiver in ("self", "cls") and cls is not None:
        return {f"{m.name}.{cls}.{attr}"} if attr in m.classes[cls] else set()
    q = symbol(receiver)
    if q is None:
        return set()
    if q in modules:                                    # module.function / module.Class
        c = _constructor(modules, f"{q}.{attr}")
        return {c} if c else set()
    module, _, name = q.rpartition(".")                 # Class.method
    methods = modules[module].classes.get(name, set()) if module in modules else set()
    return {f"{q}.{attr}"} if attr in methods else set()


def static_edges(names: tuple[str, ...] = LIVE) -> tuple[set[Edge], set[str]]:
    """All resolvable function->function edges, plus every defined function."""
    modules = {n: Module(n) for n in names}
    edges: set[Edge] = set()
    for m in modules.values():
        scopes: list[tuple[str, str | None, ast.AST]] = []
        for node in m.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                scopes.append((f"{m.name}.{node.name}", None, node))
            elif isinstance(node, ast.ClassDef):
                scopes += [(f"{m.name}.{node.name}.{f.name}", node.name, f) for f in node.body
                           if isinstance(f, (ast.FunctionDef, ast.AsyncFunctionDef))]
        for caller, cls, node in scopes:               # nested defs fold into their parent
            # "convert = TABLE[key]; convert(...)" dispatches to every table entry.
            aliases = {t.id: m.tables[n.value.value.id] for n in ast.walk(node) if isinstance(n, ast.Assign)
                       and isinstance(n.value, ast.Subscript) and isinstance(n.value.value, ast.Name)
                       and n.value.value.id in m.tables for t in n.targets if isinstance(t, ast.Name)}
            for call in ast.walk(node):
                if isinstance(call, ast.Call) and isinstance(call.func, ast.Name) and call.func.id in aliases:
                    edges |= {(caller, t) for t in aliases[call.func.id]}
                if isinstance(call, ast.Call):
                    edges |= {(caller, t) for t in _resolve(modules, m, cls, call.func) if t != caller}
    defined = set().union(*(m.functions for m in modules.values()))
    return {e for e in edges if e[1] in defined}, defined


def import_edges() -> set[Edge]:
    """Module-level dependency graph over package modules and entry points."""
    known = set(LIVE+SCRIPTS)
    sources = {n: module_path(n) for n in LIVE+SCRIPTS} | {p: ROOT/p for p in ENTRY_POINTS}
    edges: set[Edge] = set()
    for source, path in sources.items():
        package = source.rsplit(".", 1)[0] if source in known else ""
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom):
                base = (package+("."+node.module if node.module else "")) if node.level else node.module or ""
                if not node.level and package and f"{package}.{base}" in known:
                    base = f"{package}.{base}"
                targets = [base] + [f"{base}.{a.name}" for a in node.names]
            elif isinstance(node, ast.Import):
                targets = [a.name for a in node.names]
            else:
                continue
            for t in targets:
                if t in ("bcov", "quintic_bcov_lab"):
                    t += ".__init__"
                if t in known or t.endswith(".__init__"):
                    edges.add((source, t))
    # run_all_checks.py runs these scripts as subprocesses rather than importing them.
    edges |= {("run_all_checks.py", f"bcov.{s}") for s in ("test_frames", "test_gu_gauge", "bcov_yy",
              "audit_coefficients", "audit_denominators", "verify_physical_qseries")}
    # The package __init__ files re-export the core modules.
    edges |= {("bcov.__init__", "bcov.bcov_yy"), ("bcov.__init__", "bcov.quintic_frames")}
    return {e for e in edges if e[0] != e[1]}


# --------------------------------------------------------------- runtime edges

def _owner(frame, files: dict[str, str]) -> str | None:
    """Qualified name of the live function a frame belongs to.

    Comprehensions and lambdas fold into their parent; dataclass-generated
    ``__init__`` frames (filename "<string>") are skipped to reach the caller.
    """
    while frame is not None:
        filename = frame.f_code.co_filename
        if filename.startswith("<") and frame.f_code.co_name == "__init__":
            frame = frame.f_back
            continue
        module = files.get(filename)
        if module is None:
            return None
        name = frame.f_code.co_qualname
        if "<" not in name.split(".")[-1]:
            return f"{module}.{name.split('.<locals>')[0]}"
        frame = frame.f_back
    return None


@contextlib.contextmanager
def _profiled(edges: set[Edge]) -> Iterator[None]:
    files = {str(module_path(n)): n for n in LIVE}

    def profile(frame, event, arg):
        if event != "call" or "<" in frame.f_code.co_qualname:
            return
        callee = _owner(frame, files)
        caller = _owner(frame.f_back, files) if callee else None
        if caller and caller != callee:
            edges.add((caller, callee))

    sys.setprofile(profile)
    try:
        yield
    finally:
        sys.setprofile(None)


def runtime_edges() -> set[Edge]:
    """Run the test suite and the CLI_RUNS in-process under a profiler."""
    import importlib
    import pytest

    sys.path.insert(0, str(ROOT))
    edges: set[Edge] = set()
    with _profiled(edges):
        code = pytest.main(["-q", "-p", "no:cacheprovider", str(ROOT/"tests")])
    if code != 0:
        raise SystemExit(f"Test suite failed (exit {code}); call graph not regenerated")
    for module, argv in CLI_RUNS:
        main: Callable[[list[str]], int] = importlib.import_module(module).main
        with contextlib.redirect_stdout(io.StringIO()), _profiled(edges):
            main(argv)
    return edges


# ---------------------------------------------------------------------- output

def _group(name: str) -> str:
    """module.Class for a method, else the module."""
    for module in LIVE:
        if name.startswith(module+"."):
            rest = name[len(module)+1:]
            return f"{module}.{rest.split('.')[0]}" if "." in rest else module
    return name


def _label(name: str) -> str:
    """Cluster label: drop the package, keep module (and class)."""
    return name.split(".", 1)[1] if name.startswith(("bcov.", "quintic_bcov_lab.")) else name


def dot_graph(edges: dict[Edge, str], title: str, nodes: dict[str, str] | None = None) -> str:
    """Graphviz source. With nodes=None, one cluster per class/module; otherwise
    nodes maps each node to its cluster label (or "" for unclustered)."""
    names = {n for e in edges for n in e}
    nodes = nodes or {n: _label(_group(n)) for n in names}
    clusters: dict[str, list[str]] = {}
    for n in sorted(names):
        clusters.setdefault(nodes.get(n, ""), []).append(n)
    lines = [f'digraph "{title}" {{',
             '  graph [rankdir=LR, newrank=true, nodesep=0.25, ranksep=0.9, fontname="Helvetica",'
             f' label="{title}", labelloc=t, fontsize=18];',
             '  node [shape=box, style="rounded,filled", fillcolor="#ffffff", fontname="Helvetica", fontsize=11];',
             '  edge [color="#555555", arrowsize=0.6];']
    for i, (cluster, members) in enumerate(sorted(clusters.items())):
        colour = CLUSTER_COLOURS[i % len(CLUSTER_COLOURS)]
        if cluster:
            lines.append(f'  subgraph cluster_{i} {{ label="{cluster}"; style="rounded,filled"; '
                         f'fillcolor="{colour}"; color="#999999";')
        lines += [f'    "{n}" [label="{_short(n, cluster)}"];' for n in members]
        if cluster:
            lines.append("  }")
    lines += [f'  "{a}" -> "{b}" [{EDGE_STYLES[k]}];' for (a, b), k in sorted(edges.items())]
    return "\n".join(lines+["}"])


def _short(name: str, cluster: str) -> str:
    """Node label relative to its cluster; module groups say they hold module functions."""
    if name in LIVE:
        return _label(name)+" (module functions)"
    if cluster and _label(name).startswith(cluster+"."):
        return _label(name)[len(cluster)+1:]
    return _label(name)


def collapse(edges: dict[Edge, str]) -> dict[Edge, str]:
    """Class/module-level edges; a group edge is "both" if any member edge is."""
    rank = {"static": 0, "runtime": 1, "both": 2}
    out: dict[Edge, str] = {}
    for (a, b), kind in edges.items():
        key = (_group(a), _group(b))
        if key[0] != key[1] and rank[kind] >= rank[out.get(key, "static")]:
            out[key] = kind
    return out


def graphs(kinds: dict[Edge, str]) -> dict[str, tuple[str, dict[Edge, str], dict[str, str] | None]]:
    """name -> (title, edges, cluster map or None)."""
    def where(pred: Callable[[str, str], bool]) -> dict[Edge, str]:
        return {e: k for e, k in kinds.items() if pred(*e)}

    poly = lambda n: n.startswith("bcov.bcov_yy.Poly.") or n in POLY_HELPERS
    core = lambda n: n.startswith(("bcov.bcov_yy.", "bcov.quintic_frames."))
    lab = lambda n: n.startswith("quintic_bcov_lab.")
    overview = collapse(kinds)
    return {
        "overview": ("Call graph by class and module", overview,
                     {n: n.split(".")[0] for e in overview for n in e}),
        "poly": ("Poly and helpers", where(lambda a, b: poly(a) and poly(b)), None),
        "solver": ("BCOV solver (calls into Poly omitted; see overview)",
                   where(lambda a, b: a.startswith("bcov.bcov_yy.") and not poly(a) and not poly(b)), None),
        "frames": ("Frame conversions", where(lambda a, b: a.startswith("bcov.quintic_frames.")
                                               and b not in PRIMITIVE), None),
        "analysis": ("Data, patterns, physical verification and CLI (Poly primitives omitted)",
                     where(lambda a, b: not core(a) and not lab(a) and b not in PRIMITIVE
                           and not lab(b) and not (b.startswith("bcov.verify_physical_qseries.Series.")
                                                   and a.startswith("bcov.verify_physical_qseries.Series."))), None),
        "lab": ("quintic_bcov_lab research adapter", where(lambda a, b: lab(a) and b not in PRIMITIVE), None),
    }


def modules_graph() -> str:
    def cluster(n: str) -> str:
        if n.endswith(".py"):
            return "entry points (scripts, notebooks, tests)"
        if n in SCRIPTS:
            return "bcov historical scripts"
        return n.split(".")[0]
    edges = {e: "both" for e in import_edges()}          # drawn plainly: imports carry no evidence
    title = "Module import graph (live modules, historical scripts, entry points)"
    text = dot_graph(edges, title, {n: cluster(n) for e in edges for n in e})
    return text.replace("rankdir=LR", "rankdir=LR, concentrate=true")


def render(kinds: dict[Edge, str]) -> dict[str, int]:
    OUT.mkdir(parents=True, exist_ok=True)
    for stale in ("pyan3_false_positives.txt",):
        (OUT/stale).unlink(missing_ok=True)
    sources = {"modules": modules_graph()}
    for name, (title, edges, nodes) in graphs(kinds).items():
        sources[name] = dot_graph(edges, title, nodes)
    for name, text in sources.items():
        path = OUT/f"{name}.dot"
        path.write_text(text+"\n", encoding="utf-8")
        subprocess.run(["dot", "-Tsvg", str(path), "-o", str(OUT/f"{name}.svg")], check=True)
    (OUT/"edges.tsv").write_text("caller\tcallee\tevidence\n"+"".join(
        f"{a}\t{b}\t{k}\n" for (a, b), k in sorted(kinds.items())), encoding="utf-8")
    return {k: sum(v == k for v in kinds.values()) for k in EDGE_STYLES}


def update_doc(kinds: dict[Edge, str], defined: set[str], counts: dict[str, int], traced: bool) -> None:
    """Rewrite the block between the callgraph markers in ARCHITECTURE.md."""
    begin, end = "<!-- callgraph:begin -->", "<!-- callgraph:end -->"
    called = {n for e in kinds for n in e}
    rows = ["| Module | Functions/methods | Outgoing calls | Static-only (untested) |", "|---|---:|---:|---:|"]
    for module in LIVE:
        own = {f for f in defined if f.startswith(module+".")}
        out = [k for (a, _), k in kinds.items() if a in own]
        rows.append(f"| `{module}` | {len(own)} | {len(out)} | {out.count('static')} |")
    unreached = sorted(f for f in defined if f not in called)
    note = ("Runtime evidence: the full pytest suite plus " + ", ".join(
            f"`{m.split('.')[0]} {' '.join(a)}`" for m, a in CLI_RUNS) + "."
            if traced else "Generated with `--no-trace`: static edges only.")
    block = "\n".join([
        begin, "<!-- Generated by bcov/tools/callgraph.py; edits here are overwritten. -->", "",
        f"{sum(counts.values())} call edges: {counts['both']} confirmed (static and runtime), "
        f"{counts['runtime']} runtime-only, {counts['static']} static-only. {note}", "",
        *rows, "",
        f"Live functions with no call edge to or from another live function ({len(unreached)}); "
        "they are public API called only from tests, notebooks or user code:", "",
        ", ".join(f"`{_label(f)}`" for f in unreached), "",
        end])
    text = DOC.read_text(encoding="utf-8")
    head, sep, rest = text.partition(begin)
    if not sep or end not in rest:
        raise SystemExit(f"{DOC.name} lacks the {begin} ... {end} markers")
    DOC.write_text(head+block+rest.split(end, 1)[1], encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--no-trace", action="store_true", help="Skip running the tests; static edges only")
    args = parser.parse_args(argv)
    static, defined = static_edges()
    runtime = set() if args.no_trace else {e for e in runtime_edges() if e[0] in defined and e[1] in defined}
    kinds = {e: "both" for e in static & runtime}
    kinds |= {e: "runtime" for e in runtime - static}
    kinds |= {e: "static" for e in static - runtime}
    counts = render(kinds)
    update_doc(kinds, defined, counts, not args.no_trace)
    print(f"Wrote {OUT.relative_to(ROOT)}/ and updated {DOC.name}: {counts}")


if __name__ == "__main__":
    main()
