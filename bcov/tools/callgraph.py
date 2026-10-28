"""Reconciled call graph of bcov_yy and quintic_frames, rendered with Graphviz.

pyan3 alone both invents edges (name-based guesses) and misses calls made on
objects other than self. This script runs the given scripts under a profiler
to record the caller->callee edges that actually happen, merges them with the
pyan3 dot output, and writes overview/poly/solver/frames .dot and .svg files (Graphviz
`dot` must be on PATH) plus the list of rejected pyan3 edges to OUT_DIR.

    pyan3 bcov_yy.py quintic_frames.py --uses --no-defines --dot > /tmp/cg.dot
    python tools/callgraph.py /tmp/cg.dot OUT_DIR test_frames.py "bcov_yy.py --genus 4"
"""
from __future__ import annotations

import os
from pathlib import Path
import re
import runpy
import shlex
import subprocess
import sys
from collections.abc import Iterable

MODULES = {"bcov_yy.py": "bcov_yy", "quintic_frames.py": "quintic_frames"}
Edge = tuple[str, str]

# Real calls the traced scripts do not exercise. Checked by reading the code;
# the last two go through the quintic_frames._CONVERSIONS table.
UNEXERCISED: set[Edge] = {
    ("bcov_yy.Poly.__repr__", "bcov_yy.Poly.degree"),
    ("bcov_yy.Poly.evaluate", "bcov_yy.rational"),
    ("bcov_yy.Poly.evaluate", "bcov_yy.variable_index"),
    ("bcov_yy.demo", "bcov_yy.BCOV.save"),
    ("bcov_yy.BCOV.compact_recipes", "bcov_yy.RecipeArena.__init__"),
    ("bcov_yy.BCOV.compact_recipes", "bcov_yy.YYRing.clear_cache"),
    ("bcov_yy.YYRing.theta_z", "bcov_yy.YYRing.theta"),
    ("quintic_frames.FramedPolynomial.to", "quintic_frames.hkq_to_yy"),
    ("quintic_frames.FramedPolynomial.to", "quintic_frames.z_to_hkq"),
}
HELPERS = {"bcov_yy.rational", "bcov_yy.variable_index", "bcov_yy._weight", "bcov_yy._product"}
# Arithmetic/construction used everywhere; hidden as targets in the solver diagram.
PRIMITIVE = HELPERS | {f"bcov_yy.Poly.{m}" for m in (
    "__init__ __add__ __sub__ __rsub__ __neg__ __mul__ __truediv__ __pow__ __eq__ "
    "__hash__ __bool__ __len__ coerce constant scale variable terms").split()}
EDGE_STYLES = {"both": 'style=solid', "runtime": 'style=dashed, color="#1f77b4"',
               "static": 'style=dotted, color="#d62728", penwidth=1.6'}
CLUSTER_COLOURS = ("#eef3fb", "#fdf1e6", "#eaf6ec", "#f6ecf6", "#fbf8e3", "#e8f4f6", "#f3f0ea")


def _owner(frame) -> str | None:
    """Qualified name of the module function a frame belongs to (genexprs fold into their parent)."""
    while frame is not None:
        module = MODULES.get(os.path.basename(frame.f_code.co_filename))
        if module is None:
            return None
        name = frame.f_code.co_qualname
        if "<" not in name.split(".")[-1]:
            return f"{module}.{name.split('.<locals>')[0]}"
        frame = frame.f_back
    return None


def trace(scripts: Iterable[str]) -> set[Edge]:
    """Run each "script args..." as __main__ and collect caller->callee edges."""
    edges: set[Edge] = set()

    def profile(frame, event, arg):
        if event != "call" or "<" in frame.f_code.co_qualname:
            return
        callee = _owner(frame)
        caller = _owner(frame.f_back) if callee else None
        if caller and (caller != callee or frame.f_back.f_code is frame.f_code):
            edges.add((caller, callee))

    sys.path.insert(0, os.getcwd())
    for command in scripts:
        sys.argv = shlex.split(command)
        sys.setprofile(profile)
        try:
            runpy.run_path(sys.argv[0], run_name="__main__")
        except SystemExit:
            pass
        finally:
            sys.setprofile(None)
    return edges


def static_edges(dot_text: str) -> set[Edge]:
    """Function->function edges from pyan3 --dot output."""
    kinds = re.findall(r'"([^"]+)" \[label=[^\]]*?tooltip="([^\\"]+)\\n[^"]*?\\n(\w+) in ', dot_text)
    names = {node: name for node, name, kind in kinds if kind in ("function", "method")}
    return {(names[a], names[b]) for a, b in re.findall(r'"([^"]+)" -> "([^"]+)"', dot_text)
            if a in names and b in names}


def _group(name: str) -> str:
    """Class name for a method, else the module name."""
    module, rest = name.split(".", 1)
    return rest.split(".")[0] if "." in rest else module


def _short(name: str) -> str:
    """Name without module or class prefix (the cluster label shows those)."""
    return name.rsplit(".", 1)[1]


def dot_graph(edges: dict[Edge, str], title: str, clustered: bool = True) -> str:
    """Graphviz source: one labelled cluster per class/module, edge style by evidence."""
    groups: dict[str, set[str]] = {}
    for edge in edges:
        for name in edge:
            groups.setdefault(_group(name) if clustered else name, set()).add(name)
    lines = [f'digraph "{title}" {{',
             '  graph [rankdir=LR, newrank=true, nodesep=0.25, ranksep=0.9, fontname="Helvetica",'
             f' label="{title}", labelloc=t, fontsize=18];',
             '  node [shape=box, style="rounded,filled", fillcolor="#ffffff", fontname="Helvetica", fontsize=11];',
             '  edge [color="#555555", arrowsize=0.6];']
    for i, (group, names) in enumerate(sorted(groups.items())):
        if not clustered:
            lines.append(f'  "{group}" [label="{_short(group)}", fillcolor="{CLUSTER_COLOURS[i % len(CLUSTER_COLOURS)]}"];')
            continue
        lines.append(f'  subgraph cluster_{i} {{ label="{group}"; style="rounded,filled"; '
                     f'fillcolor="{CLUSTER_COLOURS[i % len(CLUSTER_COLOURS)]}"; color="#999999";')
        lines += [f'    "{n}" [label="{_short(n)}"];' for n in sorted(names)]
        lines.append("  }")
    lines += [f'  "{a}" -> "{b}" [{EDGE_STYLES[kind]}];' for (a, b), kind in sorted(edges.items())]
    lines.append("}")
    return "\n".join(lines)


def overview(edges: dict[Edge, str]) -> dict[Edge, str]:
    """Collapse to class/module level; a group edge is "both" if any member edge is."""
    rank = {"static": 0, "runtime": 1, "both": 2}
    out: dict[Edge, str] = {}
    for (a, b), kind in edges.items():
        key = (f"group.{_group(a)}", f"group.{_group(b)}")
        if key[0] != key[1] and rank[kind] >= rank[out.get(key, "static")]:
            out[key] = kind
    return out


def main(dot_path: str, out_dir: str, scripts: list[str]) -> None:
    static = static_edges(open(dot_path).read())
    runtime = trace(scripts)
    kinds = {e: "both" for e in static & runtime}
    kinds |= {e: "runtime" for e in runtime - static}
    kinds |= {e: "static" for e in UNEXERCISED}
    in_poly = lambda n: n.startswith("bcov_yy.Poly.") or n in HELPERS
    graphs = {
        "overview": ("Call graph by class and module", overview(kinds)),
        "poly": ("Poly and helpers", {e: k for e, k in kinds.items()
                                      if in_poly(e[0]) and in_poly(e[1]) and e[0] != e[1]}),
        "solver": ("BCOV solver (calls into Poly omitted; see overview)",
                   {e: k for e, k in kinds.items() if e[0].startswith("bcov_yy.")
                    and not in_poly(e[0]) and not in_poly(e[1])}),
        "frames": ("Frame conversions", {e: k for e, k in kinds.items()
                                         if e[0].startswith("quintic_frames.") and e[1] not in PRIMITIVE}),
    }
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, (title, edges) in graphs.items():
        source = out / f"{name}.dot"
        source.write_text(dot_graph(edges, title, clustered=name != "overview") + "\n", encoding="utf-8")
        subprocess.run(["dot", "-Tsvg", str(source), "-o", str(out / f"{name}.svg")], check=True)
    rejected = sorted(static - runtime - UNEXERCISED)
    (out / "pyan3_false_positives.txt").write_text(
        "".join(f"{a} -> {b}\n" for a, b in rejected), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
