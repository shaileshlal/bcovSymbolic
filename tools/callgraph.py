"""Reconciled call graph of bcov_yy and quintic_frames, rendered as Mermaid.

pyan3 alone both invents edges (name-based guesses) and misses calls made on
objects other than self. This script runs the given scripts under a profiler
to record the caller->callee edges that actually happen, merges them with the
pyan3 dot output, and writes Markdown with two Mermaid diagrams.

    pyan3 bcov_yy.py quintic_frames.py --uses --no-defines --dot > /tmp/cg.dot
    python tools/callgraph.py /tmp/cg.dot callgraph.md test_frames.py "bcov_yy.py --genus 4"
"""
from __future__ import annotations

import os
import re
import runpy
import shlex
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
    ("quintic_frames.FramedPolynomial.to", "quintic_frames.hkq_to_yy"),
    ("quintic_frames.FramedPolynomial.to", "quintic_frames.z_to_hkq"),
}
HELPERS = {"bcov_yy.rational", "bcov_yy.variable_index", "bcov_yy._weight", "bcov_yy._product"}
# Arithmetic/construction used everywhere; hidden as targets in the solver diagram.
PRIMITIVE = HELPERS | {f"bcov_yy.Poly.{m}" for m in (
    "__init__ __add__ __sub__ __rsub__ __neg__ __mul__ __truediv__ __pow__ __eq__ "
    "__hash__ __bool__ __len__ coerce constant scale variable terms").split()}
ARROWS = {"both": "-->", "runtime": "-.->", "static": "==>"}


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


def _node_id(name: str) -> str:
    return name.replace(".", "_").replace("__", "_")


def mermaid(edges: dict[Edge, str]) -> str:
    """A flowchart with one subgraph per class (or module, for free functions)."""
    groups: dict[str, set[str]] = {}
    for edge in edges:
        for name in edge:
            module, rest = name.split(".", 1)
            groups.setdefault(rest.split(".")[0] if "." in rest else module, set()).add(name)
    lines = ["```mermaid", "flowchart LR"]
    for group, names in sorted(groups.items()):
        lines.append(f"  subgraph {group}")
        lines += [f'    {_node_id(n)}["{n.split(".", 1)[1]}"]' for n in sorted(names)]
        lines.append("  end")
    lines += [f"  {_node_id(a)} {ARROWS[kind]} {_node_id(b)}" for (a, b), kind in sorted(edges.items())]
    lines.append("```")
    return "\n".join(lines)


def main(dot_path: str, out_path: str, scripts: list[str]) -> None:
    static = static_edges(open(dot_path).read())
    runtime = trace(scripts)
    kinds = {e: "both" for e in static & runtime}
    kinds |= {e: "runtime" for e in runtime - static}
    kinds |= {e: "static" for e in UNEXERCISED}
    in_poly = lambda n: n.startswith("bcov_yy.Poly.") or n in HELPERS
    poly = {e: k for e, k in kinds.items() if in_poly(e[0]) and in_poly(e[1]) and e[0] != e[1]}
    solver = {e: k for e, k in kinds.items() if not in_poly(e[0]) and e[1] not in PRIMITIVE}
    rejected = sorted(static - runtime - UNEXERCISED)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("### `Poly` and helpers\n\n" + mermaid(poly) + "\n")
        out.write("\n### Solver and frames\n\n" + mermaid(solver) + "\n")
        out.write("\npyan3 edges rejected as false positives:\n\n")
        out.write("".join(f"- `{a.split('.', 1)[1]}` → `{b.split('.', 1)[1]}`\n" for a, b in rejected))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
