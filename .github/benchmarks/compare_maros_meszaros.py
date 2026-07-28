#!/usr/bin/env python3
"""Create an aggregate current-vs-master Maros-Meszaros report."""

import argparse
import csv


TOLERANCES = ("default", "low", "med", "high")


def load(path):
    with open(path, newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    return {(row["problem"], row["tolerance"]): row for row in rows}


def percent_change(baseline, current):
    if baseline == 0.0:
        return "n/a"
    return f"{100.0 * (current - baseline) / baseline:+.1f}%"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline")
    parser.add_argument("current")
    parser.add_argument("report")
    parser.add_argument("--threshold", type=float, default=25.0)
    args = parser.parse_args()

    baseline = load(args.baseline)
    current = load(args.current)
    if baseline.keys() != current.keys():
        missing = sorted(baseline.keys() - current.keys())
        new = sorted(current.keys() - baseline.keys())
        raise RuntimeError(f"result sets differ; missing={missing}, new={new}")

    lines = [
        "## Maros–Mészáros dense benchmark",
        "",
        "Solve time is the sum of DAQP's reported solve time over every attempted "
        "problem (best of three identical runs per problem). A problem is counted "
        "as solved only when DAQP returns success and the KKT residual checks pass.",
    ]
    regression = False

    for title, predicate in (
        ("Dense subset", lambda row: True),
        ("Positive-definite subset", lambda row: row["posdef"] == "1"),
    ):
        lines.extend(
            [
                "",
                f"### {title}",
                "",
                "| Tolerance | Problems | Base solved | Current solved | Δ solved | "
                "Base solve time | Current solve time | Δ time |",
                "|:--|--:|--:|--:|--:|--:|--:|--:|",
            ]
        )
        for tolerance in TOLERANCES:
            base_rows = [
                row for (_, setting), row in baseline.items()
                if setting == tolerance and predicate(row)
            ]
            current_rows = [
                row for (_, setting), row in current.items()
                if setting == tolerance and predicate(row)
            ]
            base_solved = sum(int(row["solved"]) for row in base_rows)
            current_solved = sum(int(row["solved"]) for row in current_rows)
            base_time = sum(float(row["solve_time_s"]) for row in base_rows)
            current_time = sum(float(row["solve_time_s"]) for row in current_rows)
            time_change = 100.0 * (current_time - base_time) / base_time if base_time else 0.0
            solved_change = current_solved - base_solved
            if current_solved < base_solved or time_change > args.threshold:
                regression = True
            lines.append(
                f"| {tolerance} | {len(base_rows)} | {base_solved} | "
                f"{current_solved} | {solved_change:+d} | {base_time:.4f} s | "
                f"{current_time:.4f} s | {percent_change(base_time, current_time)} |"
            )

    lines.extend(
        [
            "",
            f"Time regressions are flagged above {args.threshold:g}%; any decrease "
            "in the solved count is also a regression.",
            "",
            "⚠️ Regression detected." if regression else "✅ No regression detected.",
        ]
    )
    report = "\n".join(lines) + "\n"
    with open(args.report, "w", encoding="utf-8") as output:
        output.write(report)
    print(report, end="")
    raise SystemExit(1 if regression else 0)


if __name__ == "__main__":
    main()
