#!/usr/bin/env python3
"""Take the screenshots of the user manual: docs/manual/img/*.png.

    python docs/manual/take_screenshots.py

Opens the GUI with the repository's data files (tle.txt, parameters.json,
ground_stations.txt, satellites.json, missions), runs a pass analysis
and the static link budgets, and grabs each tab. Needs a display (on
Linux without one: xvfb-run -s "-screen 0 1400x1000x24" python ...).
Run it again whenever the GUI layout changes, then rebuild the manual.
"""

import os
import time
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
IMG = os.path.join(HERE, "img")
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import tkinter as tk  # noqa: E402
from tkinter import ttk  # noqa: E402

from PIL import ImageGrab  # noqa: E402

import gui  # noqa: E402

WIDTH, HEIGHT = 1200, 900


def walk(widget):
    yield widget
    for child in widget.winfo_children():
        yield from walk(child)


def main() -> None:
    os.makedirs(IMG, exist_ok=True)
    gui.messagebox.showinfo = lambda *a, **k: None
    app = gui.LinkBudgetApp()
    root = app.root
    root.geometry(f"{WIDTH}x{HEIGHT}+0+0")
    root.update()

    def shot(name, height=HEIGHT):
        # Pointer on the logo (no tooltip there) and any open tooltip closed.
        root.event_generate("<Motion>", warp=True, x=40, y=40)
        root.update()
        for widget in root.winfo_children():
            if isinstance(widget, tk.Toplevel) and widget.wm_overrideredirect():
                widget.withdraw()
        root.update()
        root.after(500)
        root.update()
        x, y = root.winfo_rootx(), root.winfo_rooty()
        grab_kwargs = {"xdisplay": os.environ["DISPLAY"]} if os.environ.get("DISPLAY") else {}
        ImageGrab.grab(bbox=(x, y, x + WIDTH, y + height), **grab_kwargs).save(os.path.join(IMG, f"{name}.png"))
        print("saved", name)

    def scroll(fraction):
        target = app._wheel_targets.get(app.main_tabs.select())
        canvas = target() if target else None
        if canvas is not None:
            canvas.yview_moveto(fraction)

    def static_page(tab, text):
        for widget in walk(tab):
            if isinstance(widget, tk.Label) and widget.cget("text") == text:
                widget.event_generate("<Button-1>")
                return

    tabs = app.main_tabs.tabs()
    app.run_analysis()
    app.contact_listbox.selection_clear(0, "end")
    app.contact_listbox.selection_set(1)  # a high pass: the antenna gain changes along it
    app.on_contact_select(None)

    app.main_tabs.select(tabs[0]); scroll(0); shot("orbit_passes")
    app.main_tabs.select(tabs[1]); scroll(0); shot("downlink_params")
    scroll(1.0); shot("downlink_results")
    app.main_tabs.select(tabs[2]); app.calculate_ul_link_budget(); scroll(0); shot("uplink")

    dl = app.static_dl_tab
    app.main_tabs.select(dl)
    for page, name in (
        ("Station & Geometry", "static_geometry"),
        ("Data Channel", "static_data"),
        ("Transmitter (S/C)", "static_tx"),
        ("Propagation & PFD", "static_propagation"),
        ("Receiver (GS)", "static_rx"),
        ("Tolerances", "static_tolerances"),
    ):
        static_page(dl, page)
        shot(name)
    buttons = {w.cget("text"): w for w in walk(dl) if isinstance(w, ttk.Button)}
    buttons["Calculate"].invoke()
    shot("static_summary")
    for page, name in (
        ("Link Budget", "static_link_budget"),
        ("PFD Check", "static_pfd"),
        ("With Tolerances", "static_tolerance_results"),
    ):
        static_page(dl, page)
        shot(name)

    ul = app.static_ul_tab
    app.main_tabs.select(ul)
    static_page(ul, "Receiver (S/C)")
    shot("static_ul_rx")

    # Statistics over a few days.
    app.main_tabs.select(tabs[-1])
    app.gs_var.set("Svalbard")
    stations = app.stats_station_list
    stations.selection_clear(0, "end")
    stations.selection_set(list(stations.get(0, "end")).index("Svalbard"))
    app.stats_days_var.set("7")
    app.compute_statistics()
    while "disabled" in app.stats_button.state():
        root.update()
        time.sleep(0.1)
    root.update()
    scroll(0); shot("statistics")
    scroll(1.0); shot("statistics_tables")

    # A mission not launched yet: custom orbit and synthetic TLE.
    app.select_satellite("EPS-Sterna")
    app.main_tabs.select(tabs[0]); scroll(0); shot("custom_orbit", height=460)
    root.destroy()


if __name__ == "__main__":
    main()
