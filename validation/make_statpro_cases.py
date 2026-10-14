#!/usr/bin/env python3
"""Compute the validation cases to compare with StatPro.

Runs the Static LB · Downlink tab of the program (headless) with the AWS
parameters at a few elevations and writes validation/AWS_statpro_cases.csv:
one row per link budget line and case, with our value and an empty StatPro
column to fill in by hand. compare_statpro.py then gives the differences.

The G/T is kept constant (no sky noise model, no diameter scaling), as
StatPro has no G/T versus elevation: enter the same G/T GS in StatPro.

Usage (needs a display; on Linux without one use xvfb-run):
    python validation/make_statpro_cases.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import gui  # noqa: E402
import static_report  # noqa: E402
import tkinter as tk  # noqa: E402
from tkinter import ttk  # noqa: E402

CASES = [
    ("1", "Elevation 3°", {"elevation_deg": 3.0}),
    ("2", "Elevation 10°", {"elevation_deg": 10.0}),
    ("3", "Elevation 90° (zenith)", {"elevation_deg": 90.0}),
]
OUTPUT = os.path.join(HERE, "AWS_statpro_cases.csv")


def _walk(widget):
    yield widget
    for child in widget.winfo_children():
        yield from _walk(child)


def main() -> None:
    gui.messagebox.showinfo = lambda *a, **k: None
    gui.messagebox.showerror = lambda *a, **k: print("ERROR", *a, file=sys.stderr)
    app = gui.LinkBudgetApp()
    app.select_satellite("AWS")
    # Constant G/T, as in StatPro.
    for entry in (app.tsys_ref_entry, app.gt_ref_diameter_entry):
        entry.delete(0, tk.END)
    captured = {}
    original = app._render_fixed_link_budget_results

    def _capture(pages, budget, link_availability_pct, extra_info=None):
        captured["budget"] = budget
        captured["extra"] = extra_info or {}
        return original(pages, budget, link_availability_pct, extra_info)

    app._render_fixed_link_budget_results = _capture
    section = gui.calculations.FIXED_LINK_BUDGET_SECTIONS["downlink"]
    calculate = next(
        w for w in _walk(app.static_dl_tab) if isinstance(w, ttk.Button) and w.cget("text") == "Calculate"
    )
    base = dict(app.static_lb_values[section])
    rows = []
    for case_id, title, changes in CASES:
        values = dict(base, **changes)
        app._static_lb_appliers[section](values)
        app.root.update()
        captured.clear()
        calculate.invoke()
        app.root.update()
        budget = captured.get("budget")
        if budget is None:
            raise SystemExit(f"Case {case_id}: the static link budget was not computed (see the error above).")
        margin_req = values.get("required_margin_db")
        for row in static_report.link_budget_rows(budget, float(margin_req) if margin_req not in (None, "") else None):
            rows.append([case_id, title, row.group, row.label, f"{row.clear:.3f}", f"{row.rain:.3f}", "", ""])
        for row in static_report.pfd_rows(budget, captured["extra"]):
            rows.append([case_id, title, "PFD", row.label, f"{row.value:.3f}", "", "", ""])
        for row in static_report.tolerance_rows(budget, "rain_faded"):
            rows.append([
                case_id, title, "Tolerances (rain faded)", row.label,
                f"{row.nominal:.3f}", f"fav {row.favourable:.3f} / adv {row.adverse:.3f}", "", "",
            ])
    with open(OUTPUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Case", "Description", "Group", "Parameter",
            "Ours - Clear Sky", "Ours - Rain Faded", "StatPro - Clear Sky", "StatPro - Rain Faded",
        ])
        writer.writerows(rows)
    print(f"{len(rows)} rows written to {OUTPUT}")
    app.root.destroy()


if __name__ == "__main__":
    main()
