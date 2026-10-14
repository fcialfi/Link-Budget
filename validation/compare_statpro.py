#!/usr/bin/env python3
"""Compare our static link budget with the StatPro values written in the cases file.

Fill the StatPro columns of validation/AWS_statpro_cases.csv (Excel is
fine; keep the CSV format), then run:
    python validation/compare_statpro.py [file] [--tolerance 0.1]
Rows without a StatPro value are skipped. The script lists every
difference larger than the tolerance (dB) and a summary per case.
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _number(text):
    text = (text or "").strip().replace(",", ".")
    try:
        return float(text)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("file", nargs="?", default=os.path.join(HERE, "AWS_statpro_cases.csv"))
    parser.add_argument("--tolerance", type=float, default=0.1, help="dB (default 0.1)")
    args = parser.parse_args()
    # Excel with Italian settings saves the CSV with semicolons.
    with open(args.file, newline="", encoding="utf-8-sig") as f:
        sample = f.readline()
        f.seek(0)
        rows = list(csv.DictReader(f, delimiter=";" if sample.count(";") > sample.count(",") else ","))

    compared = differing = 0
    per_case = {}
    for row in rows:
        for column in ("Clear Sky", "Rain Faded"):
            ours = _number(row.get(f"Ours - {column}"))
            theirs = _number(row.get(f"StatPro - {column}"))
            if ours is None or theirs is None:
                continue
            compared += 1
            diff = ours - theirs
            case = per_case.setdefault(row["Case"], [0, 0, 0.0])
            case[0] += 1
            case[2] = max(case[2], abs(diff))
            if abs(diff) > args.tolerance:
                differing += 1
                case[1] += 1
                print(f"Case {row['Case']:>2}  {row['Parameter']:<45} {column:<10} "
                      f"ours {ours:9.3f}  StatPro {theirs:9.3f}  diff {diff:+.3f}")
    if not compared:
        print("No StatPro value found: fill the StatPro columns first.")
        return
    print()
    for case, (n, bad, worst) in sorted(per_case.items()):
        print(f"Case {case}: {n} values compared, {bad} differ by more than {args.tolerance} dB, "
              f"largest difference {worst:.3f} dB")
    print(f"Total: {compared} values, {differing} above the tolerance.")


if __name__ == "__main__":
    main()
