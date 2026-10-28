"""Command-line access to recursion, source patterns, and physical checks."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from importlib.resources import as_file, files
import json
from pathlib import Path

from .bcov_yy import demo
from .datasets import iter_amplitudes
from .patterns import audit_patterns
from .verify_physical_qseries import verify


def main(argv: list[str] | None = None) -> int:
    """Run a subcommand; return nonzero for failed or incomplete verification."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    recursion = commands.add_parser("recursion", help="Solve the anomaly with explicit ambiguity status")
    recursion.add_argument("--genus", type=int, default=3)
    recursion.add_argument("--checkpoint")
    patterns = commands.add_parser("patterns", help="Audit source low-X and denominator identities")
    patterns.add_argument("--genus", type=int, default=15)
    patterns.add_argument("--depth", type=int, default=4)
    patterns.add_argument("--output", type=Path)
    physical = commands.add_parser("verify", help="Compare physical q-series to GV multicover data")
    physical.add_argument("--genus", type=int, default=15)
    physical.add_argument("--degree", type=int, default=6)
    physical.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "recursion":
            demo(args.genus, args.checkpoint)
            return 0
        if not 2 <= args.genus <= 15:
            raise ValueError("Bundled amplitudes cover genera 2 through 15")
        if args.command == "patterns":
            reports = [audit_patterns(a, args.depth) for a in iter_amplitudes() if a.genus <= args.genus]
            if {r.genus for r in reports} != set(range(2, args.genus+1)):
                raise ValueError("Source genus coverage is incomplete")
            for report in reports:
                print(f"g={report.genus}: low-X={report.low_x_matches}, denominators={report.denominator_matches}")
            if args.output:
                args.output.write_text(json.dumps([asdict(r) for r in reports], default=str, indent=2)+"\n")
            return int(not all(r.low_x_matches and r.denominator_matches for r in reports))
        data = files("bcov").joinpath("data")
        with as_file(data.joinpath("bcov_hkq_atomic_genus_2_15.jsonl.gz")) as atomic:
            with as_file(data.joinpath("GVTab.csv.gz")) as gv:
                results = verify(atomic, gv, args.degree, args.genus, args.output)
        return int(any(any(x != "0" for x in r["GV_residuals"]) for r in results))
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    return 2
