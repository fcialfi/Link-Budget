"""Tkinter GUI for the satellite link budget tool."""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import tkinter.font as tkfont
import csv
import math
import os
import re
import sys
import threading
from typing import Callable

from version import APP_AUTHOR, APP_VERSION

# ---------------------------------------------------------------------------
# Visual theme
# ---------------------------------------------------------------------------
# Defined up here, before the slower imports below, because the startup
# splash screen needs it immediately.
# MATLAB-like look: light grey panels, white fields, MATLAB blue accents, a
# toolstrip-style tab bar under a dark header carrying the logo.
PALETTE = {
    "bg": "#f0f0f0",
    "surface": "#ffffff",
    "surface_alt": "#f7f7f7",
    "border": "#c8c8c8",
    "toolstrip": "#e6e6e6",
    "toolstrip_hover": "#d6d6d6",
    "header_bg": "#1f2a36",
    "header_fg": "#ffffff",
    "header_fg_muted": "#b8c4d0",
    "header_hover": "#2c3b4b",
    "accent": "#0072bd",
    "accent_dark": "#005a96",
    "accent_light": "#e3f0fa",
    "text": "#1a1a1a",
    "text_muted": "#5f5f5f",
    "danger": "#d64545",
    "success": "#2e7d32",
    "danger_dark": "#b23434",
    "row_odd": "#ffffff",
    "row_even": "#f5f9fd",
    "selection": "#cfe5f7",
    "warning": "#c77700",
    "warning_on_dark": "#ffb74d",
    "success_light": "#e8f5e9",
    "danger_light": "#fdecec",
    "danger_border": "#f1b5b5",
    "nav_disabled": "#aab3c0",
}

# MATLAB's default line colour order, for the plots.
MATLAB_COLORS = {
    "blue": "#0072BD",
    "orange": "#D95319",
    "yellow": "#EDB120",
    "purple": "#7E2F8E",
    "green": "#77AC30",
    "cyan": "#4DBEEE",
    "red": "#A2142F",
}


# Logo and icon images (assets/), found next to this file or, in the packaged
# executable, in PyInstaller's extraction folder.
_ASSETS_DIR = os.path.join(
    getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__))), "assets"
)
LOGO_HEADER_FILE = "exostaff_logo_header.png"  # 44 px high, for the window headers
LOGO_SPLASH_FILE = "exostaff_logo_splash.png"  # 60 px high, for the startup splash
LOGO_FILE = "exostaff_logo.png"  # full size, for the PDF report
ICON_FILE = "exostaff_icon.ico"
ICON_PNG_FILE = "exostaff_icon_64.png"


def _asset_path(name: str) -> str | None:
    """Path of an assets/ file, or None when it is missing."""

    path = os.path.join(_ASSETS_DIR, name)
    return path if os.path.isfile(path) else None


def _load_image(name: str, master=None) -> "tk.PhotoImage | None":
    """An assets/ PNG as a Tk image, or None when missing or unreadable.

    The caller must keep a reference to the image (e.g. on the widget that
    shows it), or Tk discards it.
    """

    path = _asset_path(name)
    if path is None:
        return None
    try:
        return tk.PhotoImage(file=path, master=master)
    except tk.TclError:
        return None


def _pick_font(preferred: list[str], fallback: str = "TkDefaultFont") -> str:
    """Return the first available font family from ``preferred``.

    Requires a Tk root to already exist; falls back silently otherwise so
    the GUI still starts on platforms missing every preferred family.
    """

    try:
        available = set(tkfont.families())
    except tk.TclError:
        return fallback
    for name in preferred:
        if name in available:
            return name
    return fallback


# ---------------------------------------------------------------------------
# Startup splash screen
# ---------------------------------------------------------------------------
# ``matplotlib``/``pandas``/``skyfield``/``astropy``/``itur`` (imported right
# below) can take several seconds to load -- especially the first time a
# frozen PyInstaller build unpacks them -- during which the app would
# otherwise show nothing at all and look hung. This splash uses only
# ``tkinter`` (already imported above) so it can appear before any of those
# heavier imports even start, then gets reused as the main window once the
# real UI is ready (see ``LinkBudgetApp._build_ui``).


def _show_startup_splash():
    """Create and immediately display a "Loading..." popup."""

    root = tk.Tk()
    root.overrideredirect(True)
    root.attributes("-topmost", True)
    root.configure(bg=PALETTE["header_bg"])

    font_family = _pick_font(["Segoe UI", "Helvetica Neue", "Helvetica", "Arial"])

    width, height = 420, 200
    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()
    root.geometry(f"{width}x{height}+{(screen_w - width) // 2}+{(screen_h - height) // 2}")

    logo = _load_image(LOGO_SPLASH_FILE, master=root)
    if logo is not None:
        logo_label = tk.Label(root, image=logo, bg=PALETTE["header_bg"])
        logo_label.image = logo
        logo_label.pack(pady=(20, 8))
    else:
        tk.Label(
            root, text="LB", bg=PALETTE["accent"], fg="white", font=(font_family, 16, "bold"), width=3
        ).pack(pady=(24, 10))
    tk.Label(
        root,
        text="Satellite Link Budget Tool",
        bg=PALETTE["header_bg"],
        fg=PALETTE["header_fg"],
        font=(font_family, 12, "bold"),
    ).pack()
    tk.Label(
        root,
        text=f"Version {APP_VERSION}  ·  Loading, please wait...",
        bg=PALETTE["header_bg"],
        fg=PALETTE["header_fg_muted"],
        font=(font_family, 9),
    ).pack(pady=(4, 2))
    tk.Label(
        root,
        text=f"Author: {APP_AUTHOR}",
        bg=PALETTE["header_bg"],
        fg=PALETTE["header_fg_muted"],
        font=(font_family, 8),
    ).pack(pady=(0, 12))

    style = ttk.Style(root)
    style.theme_use("clam")
    style.configure(
        "Splash.Horizontal.TProgressbar",
        troughcolor=PALETTE["header_bg"],
        background=PALETTE["accent"],
        bordercolor=PALETTE["header_bg"],
        lightcolor=PALETTE["accent"],
        darkcolor=PALETTE["accent"],
    )
    progress = ttk.Progressbar(
        root, mode="indeterminate", length=280, style="Splash.Horizontal.TProgressbar"
    )
    progress.pack()
    progress.start(15)

    root.update()
    return root, progress


def _splash_tick(root, progress):
    """Pump the Tk event loop once so the splash stays responsive/animated.

    Safe to call even after the splash has been torn down.
    """

    try:
        progress.step(8)
        root.update()
    except tk.TclError:
        pass


_SPLASH_ROOT, _SPLASH_PROGRESS = _show_startup_splash()

# The packaged executable shows the bootloader's splash screen from the
# moment it is started (see gui.spec); hand over to the loading window above
# now that it is up. The module only exists in a PyInstaller build with a
# splash screen.
try:
    import pyi_splash  # type: ignore[import-not-found]

    pyi_splash.close()
except Exception:
    pass

from datetime import datetime, timedelta, timezone

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MinuteLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Plots in the MATLAB style: white axes with a full box, light solid grid,
# bold titles and MATLAB's line colours.
plt.rcParams.update(
    {
        "axes.facecolor": "white",
        "axes.edgecolor": "#262626",
        "axes.linewidth": 0.8,
        "axes.titleweight": "bold",
        "axes.titlesize": 11,
        "axes.prop_cycle": plt.cycler(color=list(MATLAB_COLORS.values())),
        "grid.color": "#dcdcdc",
        "grid.linestyle": "-",
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.6,
        "legend.fontsize": 8,
        "legend.framealpha": 0.9,
    }
)

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import pandas as pd
import numpy as np

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
from skyfield.api import load, EarthSatellite, wgs84

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import astropy.units as u
import json

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import calculations
import tle_source
import missions
import licensing
import static_report

# Add path of the current script (works also in PyInstaller .exe)
if getattr(sys, 'frozen', False):
    base_path = sys._MEIPASS  # PyInstaller temp path
else:
    base_path = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, base_path)
if os.environ.get("GUI_DEBUG"):
    print(f"Running from base_path: {base_path}")
    print("Current directory content:", os.listdir(base_path))


from calculations import (
    calculate_link_budget_parameters,
    calculate_fixed_elevation_link_budget,
    atmospheric_attenuation,
    prepare_topocentric_data,
    load_antenna_pattern,
    reload_ground_stations,
    MIN_ELEVATION_DEG,
)

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)


def _format_optional(value, default: str) -> str:
    """``value`` as entry text, or ``default`` when it is missing (None)."""

    return default if value is None else f"{value:g}" if isinstance(value, (int, float)) else str(value)


def _get_optional_float(entry: ttk.Entry) -> float | None:
    """Return ``float(entry.get())`` or ``None`` when the field is empty."""

    text = entry.get().strip()
    return float(text) if text else None


def _classify_pass_direction(sat, ts, start_time: datetime, end_time: datetime) -> str:
    """Return ``A`` for ascending passes and ``D`` for descending passes."""

    t_start = ts.utc(
        start_time.year,
        start_time.month,
        start_time.day,
        start_time.hour,
        start_time.minute,
        start_time.second,
    )
    t_end = ts.utc(
        end_time.year,
        end_time.month,
        end_time.day,
        end_time.hour,
        end_time.minute,
        end_time.second,
    )
    lat_start = sat.at(t_start).subpoint().latitude.degrees
    lat_end = sat.at(t_end).subpoint().latitude.degrees
    return "A" if lat_end >= lat_start else "D"


class _Tooltip:
    """A small hover tooltip shown after a short delay over ``widget``."""

    def __init__(self, widget, text: str, delay: int = 500):
        self.widget = widget
        self.text = text
        self.delay = delay
        self._after_id: str | None = None
        self._tip: tk.Toplevel | None = None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._hide)
        widget.bind("<ButtonPress>", self._hide)

    def _schedule(self, _event=None):
        self._cancel()
        self._after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after_id is not None:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip is not None or not self.widget.winfo_viewable():
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(
            self._tip,
            text=self.text,
            background=PALETTE["header_bg"],
            foreground=PALETTE["header_fg"],
            font=("TkDefaultFont", 9),
            padx=8,
            pady=5,
            wraplength=560,
            justify="left",
        ).pack()

    def _hide(self, _event=None):
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


def _add_tooltip(widget, text: str) -> None:
    """Attach a hover tooltip with ``text`` to ``widget``."""

    _Tooltip(widget, text)


def _info_icon(parent, text: str) -> ttk.Label:
    """A small "ⓘ" label showing ``text`` when the pointer rests on it."""

    icon = ttk.Label(parent, text="ⓘ", foreground=PALETTE["accent"], cursor="question_arrow")
    # Kept on the icon so the text can be updated (e.g. with a live breakdown).
    icon.tooltip = _Tooltip(icon, text, delay=150)
    return icon


class _TreeviewRowTooltip:
    """Hover tooltip showing a per-row text for a ``ttk.Treeview``.

    ``texts`` maps item ids to the text shown while the pointer rests on
    that row; rows without an entry show nothing.
    """

    def __init__(self, tree: ttk.Treeview, texts: dict[str, str], delay: int = 400):
        self.tree = tree
        self.texts = texts
        self.delay = delay
        self._row: str | None = None
        self._after_id: str | None = None
        self._tip: tk.Toplevel | None = None
        tree.bind("<Motion>", self._on_motion, add="+")
        tree.bind("<Leave>", self._on_leave, add="+")
        tree.bind("<ButtonPress>", self._on_leave, add="+")

    def _on_motion(self, event):
        row = self.tree.identify_row(event.y)
        if row == self._row:
            return
        self._hide()
        self._row = row
        if row and self.texts.get(row):
            x, y = event.x_root + 14, event.y_root + 16
            self._after_id = self.tree.after(self.delay, lambda: self._show(row, x, y))

    def _on_leave(self, _event=None):
        self._hide()
        self._row = None

    def _show(self, row: str, x: int, y: int):
        self._after_id = None
        if self._tip is not None or not self.tree.winfo_viewable():
            return
        self._tip = tk.Toplevel(self.tree)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(
            self._tip,
            text=self.texts[row],
            background=PALETTE["header_bg"],
            foreground=PALETTE["header_fg"],
            font=("TkDefaultFont", 9),
            padx=8,
            pady=5,
            wraplength=360,
            justify="left",
        ).pack()

    def _hide(self):
        if self._after_id is not None:
            self.tree.after_cancel(self._after_id)
            self._after_id = None
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


def _set_entry_text(entry: ttk.Entry, text: str) -> None:
    """Replace an entry's text, keeping a placeholder entry's hint consistent."""

    if isinstance(entry, _PlaceholderEntry):
        entry.set_text(text)
        return
    entry.delete(0, tk.END)
    if text:
        entry.insert(0, text)


def _configure_result_tags(table: ttk.Treeview) -> None:
    """Zebra rows plus green/red text for rows checked against a requirement."""

    table.tag_configure("evenrow", background=PALETTE["row_even"])
    table.tag_configure("oddrow", background=PALETTE["row_odd"])
    table.tag_configure("pass", foreground=PALETTE["success"])
    table.tag_configure("fail", foreground=PALETTE["danger"])
    family = _pick_font(["Segoe UI", "Helvetica Neue", "Helvetica", "Arial"])
    table.tag_configure("emphasis", font=(family, 10, "bold"))
    table.tag_configure(
        "group", background=PALETTE["accent_light"], foreground=PALETTE["accent_dark"], font=(family, 9, "bold")
    )


def _format_checked(value: float, status: str | None, digits: int = 2) -> str:
    """Format a result value, with a check / cross mark when it was checked."""

    text = f"{value:.{digits}f}"
    if status == "pass":
        return f"{text}  ✓"
    if status == "fail":
        return f"{text}  ✗"
    return text


def _project_slug(project: str) -> str:
    """File name prefix from a project name ("AWS L-band" -> "AWS_L-band_"), or ""."""

    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", project.strip()).strip("_")
    return f"{slug}_" if slug else ""


def _header_logo(parent, font_family: str) -> tk.Label:
    """The company logo at the left of a navy window header ("LB" badge if missing)."""

    logo = _load_image(LOGO_HEADER_FILE, master=parent)
    if logo is not None:
        label = tk.Label(parent, image=logo, bg=PALETTE["header_bg"], bd=0)
        label.image = logo
        return label
    return tk.Label(
        parent,
        text="LB",
        bg=PALETTE["accent"],
        fg="white",
        font=(font_family, 13, "bold"),
        width=3,
        padx=4,
        pady=4,
    )


# Hover definitions for the rows of the fixed-elevation link budget results.
_FIXED_LB_ROW_HELP = {
    "EIRP [dBW]": (
        "Effective isotropic radiated power: from the Tx chain (Tx power + effective gain) when "
        "Transmitter Power is set, otherwise the EIRP of the parameters panel."
    ),
    "Free Space Loss [dB]": "Spreading loss over the slant range d: 20·log10(4π·d / λ).",
    "Gas Attenuation [dB]": "Oxygen and water-vapour absorption (ITU-R P.676) at this elevation.",
    "Cloud Attenuation [dB]": (
        "Cloud and fog attenuation (ITU-R P.840), exceeded for p% of the time. Not included in "
        "the Clear Sky column."
    ),
    "Rain Attenuation [dB]": (
        "Rain attenuation (ITU-R P.618) exceeded for p% of the time, p = 100 − link availability. "
        "Not included in the Clear Sky column."
    ),
    "Scintillation [dB]": (
        "Tropospheric scintillation fade (ITU-R P.618), included in the Rain Faded column when "
        "Scintillation is enabled."
    ),
    "Total Atmospheric Attenuation [dB]": "Combined gas, cloud, rain and scintillation attenuation (ITU-R P.618).",
    "Ionospheric Loss [dB]": (
        "Ionospheric loss: the manual value, or the scintillation fade derived from S4 "
        "(ITU-R P.531, half of the peak-to-peak fluctuation)."
    ),
    "Polarisation Mismatch Loss [dB]": (
        "Loss due to the different polarisation of Tx and Rx antennas: the entered value, or the "
        "worst case computed from the two axial ratios."
    ),
    "Multipath Losses [dB]": "Multipath allowance entered in the Propagation section.",
    "Other Attenuations [dB]": "Other static attenuation (Propagation section), subtracted from Rx Power.",
    "Total Propagation Loss [dB]": (
        "Free space loss + atmospheric attenuation + ionospheric, polarisation and multipath losses "
        "+ other attenuations."
    ),
    "Rx G/T [dB/K]": "Receiving figure of merit: antenna gain over system noise temperature.",
    "Rx Antenna Pointing Loss [dB]": "Receive antenna depointing loss: 12·(depointing angle / 3 dB beamwidth)².",
    "Rx Power [dBW]": (
        "Power at the receive antenna: EIRP − total propagation loss − multipath − other "
        "attenuations − Tx and Rx pointing losses."
    ),
    "Received S/No [dBHz]": (
        "Signal-to-noise density before implementation losses: Rx Power + G/T − k, with "
        "k = −228.6 dBW/K/Hz (Boltzmann constant)."
    ),
    "Modulation Degradation [dB]": "Modulator/filter implementation loss entered in the Receiver section.",
    "Receiver Degradation [dB]": "Demodulator implementation loss from the parameters panel.",
    "C/No after Degradations [dBHz]": "Received S/No minus modulation and receiver degradations.",
    "Transfer Frame Bit Rate [dBHz]": (
        "10·log10 of the bit rate Eb/No is referred to: the transfer frame bit rate, i.e. the "
        "coded bit rate without channel coding (= information bit rate × formatting factor)."
    ),
    "Received Eb/No [dB]": "Energy per bit over noise density: C/No after degradations − bit rate [dBHz].",
    "Required Eb/No [dB]": "Eb/No threshold of the modulation/coding for the required error rate.",
    "Data Recovery Margin [dB]": "Received Eb/No − required Eb/No.",
    "Margin – All Favourable [dB]": "Nominal margin with every tolerance at its favourable limit.",
    "Margin – All Adverse [dB]": "Nominal margin with every tolerance at its adverse limit.",
    "Margin Mean [dB]": "Nominal margin plus the mean offsets of the tolerance distributions.",
    "Margin σ [dB]": "Standard deviation of the margin: square root of the sum of the tolerance variances.",
    "MEAN – 3σ [dB]": "Statistical margin: mean − 3σ, required ≥ 0 dB (≈ 99.87% probability of closing the link).",
    "MARGIN – Worst Case RSS [dB]": (
        "Nominal margin − root-sum-square of the adverse tolerances, required ≥ 0 dB."
    ),
    "Power Flux Density [dBm/m²]": "PFD = EIRP − 10·log10(4π·d²), expressed in dBm/m².",
    "Power Flux Density [dBW/m²]": "PFD = EIRP − 10·log10(4π·d²), d = slant range, no atmospheric loss credit.",
    "Occupied Bandwidth [kHz]": "Coded bit rate × (1 + roll-off) / spectral efficiency, from the parameters panel.",
    "Power Flux Density in 4 kHz [dBW/m²/4kHz]": (
        "PFD referred to the 4 kHz reference bandwidth of the ITU limits: "
        "PFD − 10·log10(occupied bandwidth / 4 kHz)."
    ),
    "PFD Limit [dBW/m²/4kHz]": "Regulatory limit entered in the Data Rate & PFD Limit section.",
    "PFD Margin [dB]": "PFD limit − PFD in 4 kHz; must be ≥ 0 dB.",
}


class _PlaceholderEntry(ttk.Entry):
    """A ``ttk.Entry`` that shows a greyed-out example until the user types.

    ``get()`` is overridden to return ``""`` while the placeholder is
    showing, so callers can keep treating an untouched field exactly like
    an empty one -- no other code needs to know the placeholder exists.
    """

    def __init__(self, master, placeholder: str, **kwargs):
        super().__init__(master, **kwargs)
        self._placeholder = placeholder
        self._showing_placeholder = False
        self.bind("<FocusIn>", self._clear_placeholder)
        self.bind("<FocusOut>", self._show_placeholder)
        self._show_placeholder()

    def _show_placeholder(self, _event=None):
        if not ttk.Entry.get(self):
            self._showing_placeholder = True
            self.insert(0, self._placeholder)
            self.configure(foreground=PALETTE["text_muted"])

    def _clear_placeholder(self, _event=None):
        if self._showing_placeholder:
            self._showing_placeholder = False
            self.delete(0, tk.END)
            self.configure(foreground=PALETTE["text"])

    def get(self):
        return "" if self._showing_placeholder else ttk.Entry.get(self)

    def set_text(self, text: str) -> None:
        """Replace the real text; an empty ``text`` shows the placeholder again."""

        self._showing_placeholder = False
        self.delete(0, tk.END)
        self.configure(foreground=PALETTE["text"])
        if text:
            self.insert(0, text)
        else:
            try:
                focused = self.focus_get()
            except KeyError:  # focus inside a Tk-internal widget (e.g. combobox popdown)
                focused = None
            if focused is not self:
                self._show_placeholder()


def _stop_for_invalid_license(status) -> None:
    """Tell the user the licence is not valid, then exit."""

    global _SPLASH_ROOT
    root = _SPLASH_ROOT or tk.Tk()
    _SPLASH_ROOT = None
    root.withdraw()
    messagebox.showerror(
        f"{licensing.PRODUCT}: licence not valid",
        f"{status.message}\n\nContact the supplier for a new licence file.",
        parent=root,
    )
    root.destroy()
    raise SystemExit(1)


NO_PASS_SELECTED = "none (run Start/Refresh Analysis)"
# Satellite antenna pattern loaded at startup (next to the executable, in the
# working directory or next to the source), like tle.txt and parameters.json.
ANTENNA_PATTERN_FILE = "Antenna Pattern Default.txt"


class LinkBudgetApp:
    """Tkinter application for the satellite link budget tool.

    Everything that used to live in module-level globals (widgets, the
    current analysis DataFrame, contact windows, staleness flags, ...) is
    kept here as instance attributes instead, so a single object owns the
    application's state rather than a collection of free functions mutating
    shared globals.
    """

    def __init__(self):
        self.contact_windows = []
        self.df_all = pd.DataFrame()
        self.analysis_needs_refresh = True
        self.uplink_needs_refresh = False
        self.current_table_df = pd.DataFrame()
        self.current_gs_file = ""
        self.gs_menu: ttk.Combobox | None = None
        self.gs_file_var: tk.StringVar | None = None
        self.tle_file_var: tk.StringVar | None = None
        self.tle_file_path: str | None = None
        self.tle_info_label: ttk.Label | None = None
        # Satellite catalogue (satellites.json) and the satellite selected.
        self.catalogue: missions.Catalogue | None = None
        self.current_sat: missions.Satellite | None = None
        # False for a satellite without link budget parameters: contact analysis only.
        self.link_budget_enabled = True
        self.param_file_var: tk.StringVar | None = None
        self.uplink_bitrate_entry: ttk.Entry | None = None
        self.uplink_rolloff_entry: ttk.Entry | None = None
        self.uplink_overhead_entry: ttk.Entry | None = None
        self.uplink_spectral_eff_entry: ttk.Entry | None = None
        self.info_bitrate_ul_var: tk.StringVar | None = None
        self.channel_bw_ul_var: tk.StringVar | None = None
        self.uplink_table_frame: ttk.Frame | None = None
        self.uplink_recalc_button: ttk.Button | None = None
        self._mousewheel_stack: list = []
        # Inputs of the static (fixed-elevation) link budget tabs, per
        # parameters.json section, as last saved, calculated or loaded.
        self.static_lb_values: dict[str, dict] = {}
        # parameters.json section -> function refilling that static link budget tab.
        self._static_lb_appliers: dict[str, Callable[[dict], None]] = {}
        # parameters.json section -> function returning that tab's
        # elevation-independent losses, applied to the pass analysis.
        self._static_lb_loss_providers: dict[str, Callable[[], dict]] = {}
        # Functions refreshing the static tabs (e.g. after a new antenna pattern).
        self._static_lb_refreshers: list[Callable[[], None]] = []
        # str(tab frame) -> function returning the canvas the mouse wheel scrolls.
        self._wheel_targets: dict[str, Callable[[], tk.Canvas | None]] = {}
        self.current_param_file: str | None = None
        # Packaged executable only (see licensing.py): stop here with an
        # explanation if the licence is missing, invalid or expired.
        self.license_status = licensing.check_license()
        if not self.license_status.valid:
            _stop_for_invalid_license(self.license_status)

        self._build_ui()

    def run(self):
        """Enter the Tk event loop; blocks until the window is closed."""
        self.root.mainloop()

    # ------------------------------------------------------------------
    # File loading callbacks
    # ------------------------------------------------------------------

    def _apply_tle_file(self, file_path: str) -> None:
        """Parse ``file_path`` as a TLE and populate the TLE entries.

        Raises on any parsing problem; the caller decides how to surface
        that (a dialog for the manual "Load TLE from file" button, a
        silent skip for the best-effort startup auto-load).
        """

        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip()]
        if len(lines) < 2:
            raise ValueError("TLE file must contain at least two non-empty lines.")

        if lines[0].startswith("1 ") and lines[1].startswith("2 "):
            tle1, tle2 = lines[0], lines[1]
        elif len(lines) >= 3 and lines[1].startswith("1 ") and lines[2].startswith("2 "):
            tle1, tle2 = lines[1], lines[2]
        else:
            tle1, tle2 = lines[0], lines[1]
            if not (tle1.startswith("1 ") and tle2.startswith("2 ")):
                raise ValueError("TLE lines must start with '1 ' and '2 '.")

        self.tle_file_path = os.path.abspath(file_path)
        self._set_tle_lines(tle1, tle2, f"TLE: {self.tle_file_path}")

    def _set_tle_lines(self, tle1: str, tle2: str, origin: str, show_age: bool = True) -> None:
        """Fill the TLE entries and show where the TLE comes from and its age."""

        self.tle1_entry.delete(0, tk.END)
        self.tle1_entry.insert(0, tle1)
        self.tle2_entry.delete(0, tk.END)
        self.tle2_entry.insert(0, tle2)
        if self.tle_file_var is not None:
            text, old = origin, False
            try:
                if not show_age:
                    raise ValueError("synthetic TLE: no age")
                age = tle_source.age_days(tle1)
                text += f"  ·  epoch {tle_source.epoch(tle1):%Y-%m-%d %H:%M} UTC ({age:.1f} days old)"
                old = age > tle_source.MAX_AGE_DAYS
                if old:
                    text += " – old TLE, download a new one: the pass times may be off"
            except ValueError:
                pass
            self.tle_file_var.set(text)
            if self.tle_info_label is not None:
                self.tle_info_label.configure(
                    foreground=PALETTE["warning"] if old else PALETTE["text_muted"]
                )
        self.set_analysis_stale()

    def download_tle(self, auto: bool = False):
        """Download the TLE of the NORAD ID from the internet and save it in tle.txt.

        The download runs in a thread so that the window does not freeze
        while waiting for the server (up to a few tens of seconds). With
        ``auto`` (download on satellite selection) an error is shown in the
        TLE status line instead of a dialog.
        """

        text = self.norad_var.get().strip()
        try:
            norad_id = int(text)
            if norad_id <= 0:
                raise ValueError
        except ValueError:
            if not auto:
                messagebox.showerror(
                    "TLE Download", "Enter the NORAD catalogue number of the satellite (e.g. 60543)."
                )
            return
        target_path = self.tle_file_path or os.path.join(os.getcwd(), "tle.txt")

        self.tle_download_button.state(["disabled"])
        self.tle_file_var.set(f"Downloading the TLE of NORAD {norad_id}...")
        result: dict = {}

        def _work():
            try:
                result["tle"] = tle_source.fetch_tle(norad_id)
            except Exception as exc:
                result["error"] = exc

        worker = threading.Thread(target=_work, daemon=True)
        worker.start()

        def _poll():
            if worker.is_alive():
                self.root.after(150, _poll)
                return
            self.tle_download_button.state(["!disabled"])
            if self.norad_var.get().strip() != text:
                # Another satellite was selected meanwhile: just keep the file.
                if "tle" in result:
                    self._save_tle(result["tle"], target_path)
                return
            if "error" in result:
                if os.path.isfile(target_path):
                    self._apply_tle_file(target_path)
                else:
                    self._clear_tle("TLE: none loaded")
                if auto:
                    first_line = str(result["error"]).splitlines()[0]
                    self.tle_file_var.set(
                        self.tle_file_var.get() + f"  ·  download failed ({first_line.rstrip(':')})"
                    )
                    return
                messagebox.showerror(
                    "TLE Download",
                    f"{result['error']}\n\nThe TLE already loaded is kept. Without internet access, "
                    "copy the TLE from a web page (e.g. db.satnogs.org or n2yo.com) into tle.txt "
                    "and load it with \"Load TLE from file\".",
                )
                return
            tle = result["tle"]
            origin = f"TLE of {tle.name or norad_id}: downloaded from {tle.source}"
            error = self._save_tle(tle, target_path)
            if error:
                origin += f" (not saved: {error})"
            else:
                self.tle_file_path = os.path.abspath(target_path)
                origin += f", saved in {self.tle_file_path}"
            self._set_tle_lines(tle.line1, tle.line2, origin)

        self.root.after(150, _poll)

    @staticmethod
    def _save_tle(tle, path: str) -> str | None:
        """Write ``tle`` to ``path`` (creating its folder); return the error text, if any."""

        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                f.write(tle.text())
        except OSError as exc:
            return str(exc)
        return None

    def _clear_tle(self, status: str) -> None:
        """Empty the TLE entries (no TLE for the selected satellite yet)."""

        self.tle1_entry.delete(0, tk.END)
        self.tle2_entry.delete(0, tk.END)
        self.tle_file_var.set(status)
        if self.tle_info_label is not None:
            self.tle_info_label.configure(foreground=PALETTE["text_muted"])
        self.set_analysis_stale()

    # ------------------------------------------------------------------
    # Satellite selection (satellites.json and the missions folders)
    # ------------------------------------------------------------------

    def _load_catalogue(self, path: str) -> None:
        """Read the satellite catalogue and select its default satellite."""

        self.catalogue = missions.load_catalogue(path)
        labels = [sat.label for sat in self.catalogue.satellites]
        self.satellite_combo.configure(values=labels)
        self.satellite_frame_rows["select"].grid()
        default = self.catalogue.default or self.catalogue.satellites[0].name
        self.select_satellite(default)

    def _on_satellite_chosen(self, _event=None) -> None:
        index = self.satellite_combo.current()
        if self.catalogue is not None and index >= 0:
            self.select_satellite(self.catalogue.satellites[index].name)

    def select_satellite(self, name: str) -> None:
        """Make ``name`` the current satellite: its parameters, orbit and TLE.

        With ``missions/<name>/parameters.json`` the link budget parameters
        are loaded; without it only the contact analysis is available. A
        launched satellite uses its TLE (the copy kept in the mission
        folder, downloaded again when older than a day); a custom-orbit one
        a synthetic TLE built from the orbit section.
        """

        sat = self.catalogue.find(name) if self.catalogue else None
        if sat is None:
            return
        self.current_sat = sat
        self.satellite_combo.set(sat.label)
        self.norad_var.set("" if sat.norad_id is None else str(sat.norad_id))
        self.df_all = pd.DataFrame()
        self.contact_windows.clear()

        params_path = self.catalogue.parameters_path(sat)
        orbit = {}
        if params_path:
            try:
                self._apply_parameters_file(params_path)
                with open(params_path, "r", encoding="utf-8") as f:
                    orbit = json.load(f).get("orbit") or {}
                self._set_link_budget_enabled(True)
                self.mission_status_var.set(f"Link budget parameters: {params_path}")
            except Exception as exc:
                self._set_link_budget_enabled(False)
                self.mission_status_var.set(f"Parameters not valid ({exc}): contact analysis only.")
        else:
            self.current_param_file = self.catalogue.mission_file(sat, missions.PARAMETERS_FILE)
            self._set_link_budget_enabled(False)
            self.mission_status_var.set(
                f"No link budget parameters for {sat.name}: contact analysis only. To add the link budget, "
                f"put a parameters.json in {self.catalogue.mission_dir(sat)} (e.g. a copy of the AWS one)."
            )
        if not params_path or self.project_var.get().strip() == "":
            self.project_var.set(sat.name)
        # Sat-norad from the catalogue wins over the one in the parameters.
        self.norad_var.set("" if sat.norad_id is None else str(sat.norad_id))

        pattern_path = self.catalogue.mission_file(sat, missions.PATTERN_FILE)
        if not os.path.isfile(pattern_path):
            pattern_path = calculations.resolve_optional_data_file(ANTENNA_PATTERN_FILE)
        if pattern_path:
            try:
                load_antenna_pattern(pattern_path)
                self.pattern_file_var.set(f"Pattern: {os.path.abspath(pattern_path)}")
            except Exception as exc:
                self.pattern_file_var.set(f"Pattern not loaded: {exc}")

        self.tle_file_path = self.catalogue.mission_file(sat, missions.TLE_FILE)
        if sat.custom_orbit:
            self.orbit_alt_var.set(_format_optional(orbit.get("altitude_km"), "628"))
            self.orbit_ecc_var.set(_format_optional(orbit.get("eccentricity"), "0"))
            self.orbit_argp_var.set(_format_optional(orbit.get("arg_perigee_deg"), "0"))
            self.orbit_inc_var.set(_format_optional(orbit.get("inclination_deg"), ""))
            planes = orbit.get("ltan") or ["22:20"]
            planes = [str(p) for p in (planes if isinstance(planes, list) else [planes])]
            self.orbit_ltan_combo.configure(values=planes)
            self.orbit_ltan_var.set(planes[0])
            self.satellite_frame_rows["orbit"].grid()
            self.tle_download_button.state(["disabled"])
            self.generate_custom_tle(show_errors=False)
        else:
            self.satellite_frame_rows["orbit"].grid_remove()
            self.tle_download_button.state(["!disabled"])
            age = None
            if os.path.isfile(self.tle_file_path):
                try:
                    self._apply_tle_file(self.tle_file_path)
                    age = tle_source.age_days(self.tle1_entry.get().strip())
                except Exception as exc:
                    self._clear_tle(f"TLE file not valid: {exc}")
            else:
                self._clear_tle(f"No TLE of {sat.name} saved yet")
            if age is None or age > 1.0:
                self.download_tle(auto=True)
        self.set_analysis_stale()

    def _set_link_budget_enabled(self, enabled: bool) -> None:
        """Enable the link budget tabs, or leave only the contact analysis."""

        self.link_budget_enabled = enabled
        if not enabled:
            self.main_tabs.select(0)
        for tab_id in self.main_tabs.tabs()[1:]:
            self.main_tabs.tab(tab_id, state="normal" if enabled else "disabled")
        if enabled:
            self.contact_mode_var.set("")
        else:
            self.contact_mode_var.set(
                "Contact analysis only: the passes are listed with their duration and maximum elevation."
            )

    def _custom_orbit(self) -> dict:
        """The custom orbit inputs; raise ValueError when one is not valid."""

        altitude = float(self.orbit_alt_var.get())
        ltan = tle_source.ltan_hours(self.orbit_ltan_var.get())
        inclination_text = self.orbit_inc_var.get().strip()
        inclination = float(inclination_text) if inclination_text else None
        eccentricity = float(self.orbit_ecc_var.get().strip() or 0.0)
        arg_perigee = float(self.orbit_argp_var.get().strip() or 0.0)
        return {
            "altitude_km": altitude,
            "ltan_h": ltan,
            "inclination_deg": inclination,
            "eccentricity": eccentricity,
            "arg_perigee_deg": arg_perigee,
        }

    def generate_custom_tle(self, show_errors: bool = True) -> bool:
        """Build the synthetic TLE of the custom orbit, with epoch at the analysed date."""

        try:
            orbit = self._custom_orbit()
            try:
                day = datetime.strptime(self.date_entry.get().strip(), "%Y-%m-%d")
            except ValueError:
                day = datetime.now(timezone.utc)
            when = datetime(day.year, day.month, day.day, tzinfo=timezone.utc)
            name = self.current_sat.name if self.current_sat else "CUSTOM ORBIT"
            tle = tle_source.synthetic_tle(
                orbit["altitude_km"],
                orbit["ltan_h"],
                when,
                orbit["inclination_deg"],
                name=name.upper(),
                eccentricity=orbit["eccentricity"],
                arg_perigee_deg=orbit["arg_perigee_deg"],
            )
        except ValueError as exc:
            if show_errors:
                messagebox.showerror("Custom Orbit", f"Invalid orbit: {exc}")
            self._clear_tle(f"Custom orbit not valid: {exc}")
            return False
        inclination = float(tle.line2[8:16])
        kind = "sun-synchronous" if orbit["inclination_deg"] is None else "inclination given"
        # Kept in the mission folder, like the TLE downloaded for a launched satellite.
        tle.name = f"{tle.name} (synthetic, LTAN {self.orbit_ltan_var.get().strip()})"
        saved = ""
        if self.tle_file_path:
            error = self._save_tle(tle, self.tle_file_path)
            saved = f" (not saved: {error})" if error else f", saved in {self.tle_file_path}"
        self._set_tle_lines(
            tle.line1,
            tle.line2,
            f"Synthetic TLE{saved}: mean altitude {orbit['altitude_km']:g} km, e = {orbit['eccentricity']:g}, "
            f"i = {inclination:.3f}° ({kind}), LTAN {self.orbit_ltan_var.get().strip()}, "
            f"epoch {when:%Y-%m-%d}, nodal period {1440.0 / float(tle.line2[52:63]):.3f} min",
            show_age=False,
        )
        return True

    def _custom_orbit_for_file(self) -> dict | None:
        """The custom orbit as saved in parameters.json, or None for a launched satellite."""

        if not (self.current_sat and self.current_sat.custom_orbit):
            return None
        try:
            orbit = self._custom_orbit()
        except ValueError:
            return None
        ltan = self.orbit_ltan_var.get().strip()
        planes = [ltan] + [p for p in self.orbit_ltan_combo.cget("values") if p != ltan]
        return {
            "altitude_km": orbit["altitude_km"],
            "eccentricity": orbit["eccentricity"],
            "arg_perigee_deg": orbit["arg_perigee_deg"],
            "inclination_deg": orbit["inclination_deg"],
            "ltan": planes,
        }

    def load_tle_from_file(self):
        """Let the user pick a TLE file and populate the GUI fields."""

        file_path = filedialog.askopenfilename(
            title="Select TLE file",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
        )
        if not file_path:
            return

        try:
            self._apply_tle_file(file_path)
        except Exception as exc:
            messagebox.showerror("TLE Error", f"Failed to load TLE: {exc}")

    def _apply_parameters_file(self, file_path: str) -> None:
        """Parse ``file_path`` as a parameters JSON and fill matching entries.

        The JSON structure should provide keys such as ``eirp_sat_dbw`` or
        ``frequency_ghz``. See the README for the full schema and example.
        Missing keys are ignored so the user can provide only the fields
        they need. The optional ``fixed_elevation_downlink`` /
        ``fixed_elevation_uplink`` sections hold the inputs of the static
        (fixed-elevation) link budget tabs, which are refilled from them.
        Raises on any parsing problem, mirroring :meth:`_apply_tle_file`.
        """

        with open(file_path, "r", encoding="utf-8") as f:
            payload = json.load(f)

        field_map = self._parameter_entries()

        for key, entry in field_map.items():
            if key in payload and entry is not None:
                entry.delete(0, tk.END)
                # null leaves the field empty (e.g. C/Io: no interferer).
                if payload[key] is not None:
                    entry.insert(0, str(payload[key]))

        if "project_name" in payload:
            self.project_var.set(str(payload["project_name"] or ""))
        if "norad_id" in payload:
            self.norad_var.set("" if payload["norad_id"] is None else str(payload["norad_id"]))

        for section in calculations.FIXED_LINK_BUDGET_SECTIONS.values():
            if isinstance(payload.get(section), dict):
                self.static_lb_values[section] = payload[section]
                if section in self._static_lb_appliers:
                    self._static_lb_appliers[section](payload[section])
        self.current_param_file = os.path.abspath(file_path)

        if self.param_file_var is not None:
            self.param_file_var.set(f"Parameters: {os.path.abspath(file_path)}")

        self.update_link_budget_derived()
        self.set_analysis_stale()

    def _parameter_entries(self) -> dict:
        """``parameters.json`` key -> parameters panel entry."""

        return {
            "eirp_sat_dbw": self.eirp_sat_entry,
            "eirp_gs_dbw": self.eirp_gs_entry,
            "frequency_ghz": self.freq_entry,
            "uplink_frequency_ghz": self.uplink_freq_entry,
            "c_io_dbhz": self.CIo_entry,
            "gt_gs_dbk": self.gt_gs_entry,
            "gt_sat_dbk": self.gt_sat_entry,
            "antenna_diameter_m": self.d_gs_entry,
            "link_availability_pct": self.LA_entry,
            "other_attenuations_db": self.other_att_entry,
            "uplink_other_attenuations_db": self.other_att_ul_entry,
            "bitrate_mbps": self.bitrate_entry,
            "rolloff": self.rolloff_entry,
            "demod_loss_db": self.demod_loss_entry,
            "demod_loss_sat_db": self.demod_loss_ul_entry,
            "overhead": self.overhead_entry,
            "spectral_efficiency_bpshz": self.spectral_eff_entry,
            "uplink_bitrate_mbps": self.uplink_bitrate_entry,
            "uplink_rolloff": self.uplink_rolloff_entry,
            "uplink_overhead": self.uplink_overhead_entry,
            "uplink_spectral_efficiency_bpshz": self.uplink_spectral_eff_entry,
        }

    def _linked_var(self, entry: ttk.Entry) -> tk.StringVar:
        """A StringVar shared by a parameters panel entry and its copies.

        Created on first use and attached to ``entry`` as its textvariable,
        so another widget using the same variable edits the panel value
        directly, and the other way round.
        """

        var = getattr(entry, "_linked_var", None)
        if var is None:
            var = tk.StringVar(value=entry.get())
            entry.configure(textvariable=var)
            entry._linked_var = var
        return var

    def load_parameters_from_file(self):
        """Let the user pick a parameters JSON file and load it."""

        file_path = filedialog.askopenfilename(
            title="Select Parameters File",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
        )
        if not file_path:
            return

        try:
            self._apply_parameters_file(file_path)
        except Exception as exc:
            messagebox.showerror("Parameters", f"Unable to load parameters: {exc}")

    def load_ground_stations_from_file(self):
        """Let the user pick a ground station catalogue and refresh the UI."""

        file_path = filedialog.askopenfilename(
            title="Select Ground Stations File",
            filetypes=[("Text Files", "*.txt"), ("CSV Files", "*.csv"), ("All Files", "*.*")],
        )
        if not file_path:
            return

        try:
            stations = reload_ground_stations(file_path)
        except Exception as exc:
            messagebox.showerror("Ground Stations", f"Unable to load ground stations: {exc}")
            return

        self.current_gs_file = os.path.abspath(file_path)
        if self.gs_file_var is not None:
            self.gs_file_var.set(f"Ground stations: {self.current_gs_file}")
        if self.gs_menu is not None:
            self.gs_menu["values"] = list(stations.keys())
        if stations:
            selected = self.gs_var.get()
            if selected not in stations:
                self.gs_var.set(next(iter(stations)))
        self.set_analysis_stale()

    def _load_startup_defaults(self):
        """Auto-load a default TLE and parameters file at startup.

        Mirrors how the ground station catalogue is already picked up
        automatically from ``ground_stations.txt``: on launch, look for a
        ``tle.txt`` and a ``parameters.json`` next to the executable / in
        the working directory / next to the source, and pre-fill the TLE
        and parameter fields if found. Either file is entirely optional --
        a missing or invalid one is silently skipped so it never blocks
        startup, and the user can still override everything from the
        "Load TLE from file" / "Import Link Budget Parameters" buttons.
        """

        # The ground station catalogue is read when calculations is imported.
        if self.gs_file_var is not None:
            gs_path = calculations.GROUND_STATIONS_FILE
            self.gs_file_var.set(
                f"Ground stations: {os.path.abspath(gs_path)}"
                if gs_path and os.path.isfile(gs_path)
                else "Ground stations: built-in list (ground_stations.txt not found)"
            )

        catalogue_path = calculations.resolve_optional_data_file(
            missions.CATALOGUE_FILE, env_var="SATELLITES_FILE"
        )
        tle_path = calculations.resolve_optional_data_file("tle.txt", env_var="TLE_FILE")
        # (The satellite antenna pattern file is picked up the same way.)
        if tle_path and not catalogue_path:
            try:
                self._apply_tle_file(tle_path)
            except Exception as exc:
                print(f"Warning: could not auto-load default TLE '{tle_path}': {exc}", file=sys.stderr)

        pattern_path = calculations.resolve_optional_data_file(
            ANTENNA_PATTERN_FILE, env_var="ANTENNA_PATTERN_FILE"
        )
        if pattern_path:
            try:
                load_antenna_pattern(pattern_path)
                self.pattern_file_var.set(f"Pattern: {os.path.abspath(pattern_path)}")
            except Exception as exc:
                print(
                    f"Warning: could not auto-load the antenna pattern '{pattern_path}': {exc}",
                    file=sys.stderr,
                )

        if catalogue_path:
            try:
                self._load_catalogue(catalogue_path)
                return
            except Exception as exc:
                print(f"Warning: could not load the satellite list '{catalogue_path}': {exc}", file=sys.stderr)

        # Without satellites.json: a single tle.txt and parameters.json, as before.
        params_path = calculations.resolve_optional_data_file(
            "parameters.json", env_var="PARAMETERS_FILE"
        )
        if params_path:
            try:
                self._apply_parameters_file(params_path)
            except Exception as exc:
                print(
                    f"Warning: could not auto-load default parameters '{params_path}': {exc}",
                    file=sys.stderr,
                )

    def _bind_scroll_wheel(
        self,
        canvas: tk.Canvas,
        owner: tk.Widget | None = None,
        target: Callable[[], tk.Canvas | None] | None = None,
    ) -> None:
        """Route mouse-wheel scrolling to ``canvas``.

        Mouse-wheel events are delivered to a single widget determined by
        the platform (not simply "whatever the pointer is over"), so a
        scrollable canvas has to grab them globally via ``bind_all``. That
        would let a second scrollable area (e.g. a popup opened on top of
        the main window) permanently steal wheel-scrolling from the first
        one -- an ``<Enter>``/``<Leave>`` toggle doesn't fix this either,
        since those fire on every child-widget boundary, not just the
        canvas's own edge. Instead this keeps a small stack: binding pushes
        this canvas on top, and if ``owner`` is given (a popup
        ``Toplevel``), its ``<Destroy>`` pops this canvas back off so
        wheel-scrolling reverts to whichever one -- the main window or an
        outer popup -- was active underneath it.

        ``target``, if given, returns the canvas to scroll at the time of
        the event (e.g. the one on the page currently shown), instead of
        ``canvas``.
        """

        wheel_event = "<Mousewheel>" if sys.platform == "darwin" else "<MouseWheel>"

        def _scroll(event):
            scrolled = target() if target is not None else canvas
            if scrolled is None:
                return
            if sys.platform == "darwin":
                scrolled.yview_scroll(-1 * event.delta, "units")
            else:
                scrolled.yview_scroll(int(-1 * (event.delta / 120)), "units")

        self._mousewheel_stack.append(_scroll)
        canvas.bind_all(wheel_event, _scroll)

        if owner is not None:
            def _on_destroy(event):
                if event.widget is not owner:
                    return
                if _scroll in self._mousewheel_stack:
                    self._mousewheel_stack.remove(_scroll)
                if self._mousewheel_stack:
                    canvas.bind_all(wheel_event, self._mousewheel_stack[-1])
                else:
                    canvas.unbind_all(wheel_event)

            owner.bind("<Destroy>", _on_destroy)

    def _make_scrollable_window_body(self, win: tk.Toplevel) -> ttk.Frame:
        """Wrap ``win``'s content area in a vertically scrollable canvas.

        Returns the frame that should be used as the parent for the
        window's widgets -- everything packed/gridded into it becomes
        scrollable, with a vertical scrollbar plus mouse-wheel support
        (see :meth:`_bind_scroll_wheel`).
        """

        canvas_frame, body, canvas = self._make_scrollable_frame(win)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        self._bind_scroll_wheel(canvas, owner=win)
        return body

    def _make_scrollable_frame(self, parent: tk.Widget) -> tuple[ttk.Frame, ttk.Frame, tk.Canvas]:
        """A vertically scrollable area inside ``parent``.

        Returns ``(container, body, canvas)``: the caller places
        ``container`` and puts its widgets in ``body``; mouse-wheel
        scrolling is left to the caller (see :meth:`_bind_scroll_wheel`).
        The scrollbar is shown only when the content is taller than the
        area.
        """

        container = ttk.Frame(parent)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        canvas = tk.Canvas(container, bg=PALETTE["bg"], highlightthickness=0)
        canvas.grid(row=0, column=0, sticky="nsew")
        vscroll = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        vscroll.grid(row=0, column=1, sticky="ns")

        def _set_scrollbar(first, last):
            vscroll.set(first, last)
            if float(first) <= 0.0 and float(last) >= 1.0:
                vscroll.grid_remove()
            else:
                vscroll.grid()

        canvas.configure(yscrollcommand=_set_scrollbar)

        body = ttk.Frame(canvas)
        window_id = canvas.create_window((0, 0), window=body, anchor="nw")

        def _on_canvas_configure(event):
            canvas.itemconfig(window_id, width=event.width)
            canvas.configure(scrollregion=canvas.bbox("all"))

        canvas.bind("<Configure>", _on_canvas_configure)
        body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        return container, body, canvas

    # ------------------------------------------------------------------
    # Staleness tracking
    # ------------------------------------------------------------------

    def clear_uplink_table(self):
        """Remove any existing uplink-only table widgets."""

        if self.uplink_table_frame is None:
            return
        for widget in self.uplink_table_frame.winfo_children():
            widget.destroy()

    def set_analysis_stale(self):
        """Mark the analysis as stale so the user must recompute."""
        if not self.analysis_needs_refresh:
            self.analysis_needs_refresh = True
            self.start_refresh_button.config(style="Accent.TButton")
            self.recalc_button.config(state="disabled")
            self.clear_plot_and_table()
            self.clear_uplink_table()
            self.contact_listbox.delete(0, tk.END)
            self.selected_pass_var.set(NO_PASS_SELECTED)

    def set_uplink_budget_stale(self):
        """Mark the uplink-only link budget as stale."""

        if not self.uplink_needs_refresh:
            self.uplink_needs_refresh = True
            self.clear_uplink_table()
            if self.uplink_recalc_button is not None:
                self.uplink_recalc_button.config(style="Accent.TButton")

    def set_uplink_budget_fresh(self):
        """Reset the uplink-only link budget stale flag."""

        if self.uplink_needs_refresh:
            self.uplink_needs_refresh = False
            if self.uplink_recalc_button is not None:
                self.uplink_recalc_button.config(style="TButton")

    def _pass_loss_settings(self) -> dict:
        """Values the pass analysis takes from the static link budget tabs,
        as keyword arguments of :func:`calculate_link_budget_parameters`.

        The satellite antenna pattern is always applied: on the downlink the
        EIRP is the Tx chain of the Static LB · Downlink tab plus the pattern
        gain towards the station (or EIRP SAT at the pattern peak when that
        tab uses an EIRP value); on the uplink G/T SAT is taken at the
        pattern peak. The ground station antenna tracks the satellite, so
        it has no pointing loss. Also writes what is applied into the
        Downlink and Uplink tabs.
        """

        kwargs: dict = {"apply_antenna_pattern": True}
        directions = (
            ("Downlink", "downlink", self.include_static_dl_var, self.dl_losses_var, ""),
            ("Uplink", "uplink", self.include_static_ul_var, self.ul_losses_var, "ul_"),
        )
        for title, key, include_var, text_var, prefix in directions:
            if title == "Downlink":
                antenna = (
                    "EIRP: EIRP SAT at the peak of the satellite antenna pattern, minus the pattern "
                    "roll-off towards the station."
                )
            else:
                antenna = (
                    "G/T SAT at the peak of the satellite antenna pattern, minus the pattern roll-off "
                    "towards the station."
                )
            if not include_var.get():
                text_var.set(
                    f"{antenna}\nStatic link budget losses not included: only the losses of this tab "
                    "(atmospheric, other attenuations, demodulator loss); Eb/No referred to the "
                    "information bit rate."
                )
                continue
            provider = self._static_lb_loss_providers.get(calculations.FIXED_LINK_BUDGET_SECTIONS[key])
            try:
                values = provider() if provider else None
            except ValueError as exc:
                text_var.set(f"{antenna}\nNot applied: an input of the Static LB · {title} tab is not valid ({exc})")
                continue
            if values is None:
                continue
            if title == "Downlink" and values["tx_chain_dbw"] is not None:
                kwargs["tx_chain_dbw"] = values["tx_chain_dbw"]
                antenna = (
                    f"EIRP: Tx chain of the Static LB · Downlink tab ({values['tx_chain_dbw']:.2f} dBW: "
                    "power − circuit − VSWR losses) + satellite antenna gain from the pattern at each "
                    "time step."
                )
            parts = [
                ("ionospheric", values["ionospheric_loss_db"]),
                ("polarisation", values["polarisation_loss_db"]),
                ("multipath", values["multipath_loss_db"]),
            ]
            total = sum(value for _name, value in parts)
            kwargs[f"{prefix}link_losses_db"] = total
            kwargs[f"{prefix}modulation_degradation_db"] = values["modulation_degradation_db"]
            kwargs[f"{prefix}formatting_overhead"] = values["formatting_overhead"]
            text_var.set(
                f"{antenna}\n"
                + "Link losses: "
                + " + ".join(f"{name} {value:.2f}" for name, value in parts)
                + f" = {total:.2f} dB (GS antenna tracking: no pointing loss).   "
                + f"Modulation degradation: {values['modulation_degradation_db']:.2f} dB.   "
                + "Eb/No referred to the transfer frame bit rate (formatting factor "
                + f"{values['formatting_overhead']:.6g})."
            )
        return kwargs

    def _on_pass_losses_changed(self) -> None:
        """A loss checkbox changed: show what is applied and flag the budgets."""

        self._pass_loss_settings()
        self.set_uplink_budget_stale()

    # ------------------------------------------------------------------
    # Main Analysis Function
    # ------------------------------------------------------------------

    def run_analysis(self):
        """Run a full orbital and link budget analysis using GUI inputs.

        The function reads TLE data, ground station parameters and link
        budget settings from the GUI fields. It then propagates the orbit
        for a 24 hour period with a 30 second step, computes the link
        budget for each step using :func:`calculate_link_budget_parameters`
        and stores the results in :attr:`df_all`. Detected contact windows
        are populated in :attr:`contact_windows` and displayed in the GUI.

        No parameters are accepted; all values are taken from the GUI
        widgets. The method updates instance state and returns ``None``.
        """

        if self.current_sat is not None and self.current_sat.custom_orbit:
            # The synthetic TLE has its epoch at the analysed day.
            if not self.generate_custom_tle():
                return
        tle1 = self.tle1_entry.get().strip()
        tle2 = self.tle2_entry.get().strip()
        gs_name = self.gs_var.get()
        date_str = self.date_entry.get().strip()

        if not all([tle1, tle2, gs_name, date_str]):
            messagebox.showerror("Error", "All TLE, Ground Station, and Date fields must be filled.")
            return

        try:
            obs_date = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        except ValueError:
            messagebox.showerror("Error", "Invalid date format. Please use YYYY-MM-DD.")
            return

        if not self.link_budget_enabled:
            self._run_contact_analysis(tle1, tle2, gs_name, obs_date)
            return

        try:
            freq = float(self.freq_entry.get()) * u.GHz
            link_availability = float(self.LA_entry.get())
            p = 100.0 - link_availability
            d_gs = float(self.d_gs_entry.get())
            other_att = float(self.other_att_entry.get() or 0.0)
            eirp_sat = float(self.eirp_sat_entry.get())
            gt_gs = float(self.gt_gs_entry.get())
            cisat_db = float(self.CIo_entry.get()) if self.CIo_entry.get().strip() else None
            cisat_lin = 10 ** (cisat_db / 10.0) if cisat_db is not None else None
            bitrate = float(self.bitrate_entry.get()) * 1e6
            demod_loss = float(self.demod_loss_entry.get())
            overhead = float(self.overhead_entry.get())
            uplink_freq_val = _get_optional_float(self.uplink_freq_entry)
            uplink_freq = uplink_freq_val * u.GHz if uplink_freq_val is not None else None
            eirp_gs = _get_optional_float(self.eirp_gs_entry)
            gt_sat = _get_optional_float(self.gt_sat_entry)
            demod_loss_ul = _get_optional_float(self.demod_loss_ul_entry)
            other_att_ul = _get_optional_float(self.other_att_ul_entry)
            uplink_bitrate_val = (
                _get_optional_float(self.uplink_bitrate_entry) if self.uplink_bitrate_entry else None
            )
            uplink_overhead = (
                _get_optional_float(self.uplink_overhead_entry) if self.uplink_overhead_entry else None
            )
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid numerical input: {e}.")
            return

        lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]
        loss_kwargs = self._pass_loss_settings()
        alt_gs_km = alt_gs_m / 1000.0

        ts = load.timescale()
        try:
            sat = EarthSatellite(tle1, tle2, "SAT", ts)
        except Exception as e:
            messagebox.showerror("TLE Error", f"Invalid TLE data: {e}.")
            return
        gs = wgs84.latlon(lat_gs, lon_gs, elevation_m=alt_gs_m)

        start_time = datetime(obs_date.year, obs_date.month, obs_date.day, 0, 0, 0, tzinfo=timezone.utc)
        times = [start_time + timedelta(seconds=30 * i) for i in range(2880)]
        sky_times = ts.utc(
            [t.year for t in times],
            [t.month for t in times],
            [t.day for t in times],
            [t.hour for t in times],
            [t.minute for t in times],
            [t.second for t in times],
        )

        altitudes, slant_ranges, dopplers, off_boresight = prepare_topocentric_data(
            sat, gs, sky_times, freq
        )
        atm_att = atmospheric_attenuation(
            lat_gs,
            lon_gs,
            freq.to(u.GHz).value,
            p,
            d_gs,
            alt_gs_km,
            include_scintillation=self.scint_var.get(),
        )
        uplink_atm_att = None
        if uplink_freq is not None and eirp_gs is not None and gt_sat is not None:
            uplink_atm_att = atmospheric_attenuation(
                lat_gs,
                lon_gs,
                uplink_freq.to(u.GHz).value,
                p,
                d_gs,
                alt_gs_km,
                include_scintillation=self.scint_var.get(),
            )
        self.atm_label_var.set(
            f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}° El: {atm_att:.3f}"
        )
        results_list = []
        self.contact_windows.clear()
        in_contact = False
        contact_start_time = None

        for i, t_utc_dt in enumerate(times):
            elev = altitudes[i]
            t_sky = sky_times[i]
            params = calculate_link_budget_parameters(
                t_sky,
                sat,
                gs,
                freq,
                p,
                d_gs,
                alt_gs_km,
                eirp_sat,
                gt_gs,
                demod_loss,
                bitrate,
                overhead,
                cisat_lin,
                other_att,
                atm_att,
                elev,
                slant_ranges[i],
                dopplers[i],
                off_boresight[i],
                uplink_freq,
                eirp_gs,
                gt_sat,
                demod_loss_ul,
                other_att_ul,
                uplink_atm_att,
                uplink_bitrate_val * 1e6 if uplink_bitrate_val is not None else None,
                uplink_overhead,
                **loss_kwargs,
            )
            rounding_rules = {"Atmospheric Att (dB)": 3, "UL Atmospheric Att (dB)": 3}
            for key in [
                "Elevation (°)",
                "Slant Range (km)",
                "Path Loss (dB)",
                "Atmospheric Att (dB)",
                "S/C Antenna Gain (dBi)",
                "EIRP (dBW)",
                "S/C Antenna Loss (dB)",
                "Off Boresight Angle (°)",
                "Link Losses (dB)",
                "UL Link Losses (dB)",
                "Rx Power (dBW)",
                "C/(No+Io) (dBHz)",
                "Eb/No (dB)",
                "Doppler Shift (kHz)",
                "UL Path Loss (dB)",
                "UL Atmospheric Att (dB)",
                "UL S/C Antenna Loss (dB)",
                "UL Rx Power (dBW)",
                "UL C/No (dBHz)",
                "UL Eb/No (dB)",
            ]:
                if params.get(key) is not None:
                    decimals = rounding_rules.get(key, 2)
                    params[key] = round(params[key], decimals)
            results_list.append(params)

            is_visible_now = elev >= MIN_ELEVATION_DEG
            if is_visible_now and not in_contact:
                contact_start_time = t_utc_dt
                in_contact = True
            elif not is_visible_now and in_contact:
                contact_end_time = t_utc_dt
                pass_direction = _classify_pass_direction(sat, ts, contact_start_time, contact_end_time)
                self.contact_windows.append((contact_start_time, contact_end_time, pass_direction))
                in_contact = False
        if in_contact:
            pass_direction = _classify_pass_direction(sat, ts, contact_start_time, times[-1])
            self.contact_windows.append((contact_start_time, times[-1], pass_direction))

        self.df_all = pd.DataFrame(results_list)
        self.display_contacts(self.contact_windows)
        if self.contact_windows:
            self.contact_listbox.selection_set(0)
            self.on_contact_select(None)
        else:
            self.clear_plot_and_table()
            self.selected_pass_var.set("no pass above the minimum elevation on this date")

        self.analysis_needs_refresh = False
        self.start_refresh_button.config(style="TButton")
        self.recalc_button.config(state="!disabled")

    def _run_contact_analysis(self, tle1: str, tle2: str, gs_name: str, obs_date) -> None:
        """Passes of a satellite without link budget parameters: geometry only."""

        lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]
        ts = load.timescale()
        try:
            sat = EarthSatellite(tle1, tle2, "SAT", ts)
        except Exception as e:
            messagebox.showerror("TLE Error", f"Invalid TLE data: {e}.")
            return
        gs = wgs84.latlon(lat_gs, lon_gs, elevation_m=alt_gs_m)
        start_time = datetime(obs_date.year, obs_date.month, obs_date.day, tzinfo=timezone.utc)
        times = [start_time + timedelta(seconds=30 * i) for i in range(2880)]
        sky_times = ts.from_datetimes(times)
        altitudes, slant_ranges, _, _ = prepare_topocentric_data(sat, gs, sky_times, 1.0 * u.GHz)
        self.df_all = pd.DataFrame(
            {
                "Time (UTC)": times,
                "Elevation (°)": np.round(altitudes, 2),
                "Slant Range (km)": np.round(slant_ranges, 2),
            }
        )
        self.contact_windows.clear()
        visible = altitudes >= MIN_ELEVATION_DEG
        start = None
        for i, now_visible in enumerate(visible):
            if now_visible and start is None:
                start = times[i]
            elif not now_visible and start is not None:
                self.contact_windows.append((start, times[i], _classify_pass_direction(sat, ts, start, times[i])))
                start = None
        if start is not None:
            self.contact_windows.append((start, times[-1], _classify_pass_direction(sat, ts, start, times[-1])))
        self.atm_label_var.set("")
        self.display_contacts(self.contact_windows)
        self.selected_pass_var.set(
            f"{len(self.contact_windows)} passes above {MIN_ELEVATION_DEG:g}° (contact analysis only)"
            if self.contact_windows
            else "no pass above the minimum elevation on this date"
        )
        self.analysis_needs_refresh = False
        self.start_refresh_button.config(style="TButton")

    # ------------------------------------------------------------------
    # Recalculate Link Budget Function
    # ------------------------------------------------------------------

    def recalculate_link_budget(self):
        """Recompute the link budget for previously generated contact times.

        This method does not re-run the orbital propagation. Instead it
        uses the existing time stamps stored in :attr:`df_all` and
        recalculates the link budget with the current parameters from the
        GUI. The updated results replace the contents of :attr:`df_all` and
        the plots/tables shown in the interface are refreshed accordingly.

        All configuration values are taken from the GUI widgets and no
        value is returned.
        """
        if not self.link_budget_enabled:
            return
        if self.df_all.empty:
            messagebox.showwarning("Warning", "Please generate passes first.")
            return
        if self.analysis_needs_refresh:
            messagebox.showwarning("Warning", "Please refresh the analysis first.")
            return
        try:
            freq = float(self.freq_entry.get()) * u.GHz
            link_availability = float(self.LA_entry.get())
            p = 100.0 - link_availability
            d_gs = float(self.d_gs_entry.get())
            other_att = float(self.other_att_entry.get() or 0.0)
            eirp_sat = float(self.eirp_sat_entry.get())
            gt_gs = float(self.gt_gs_entry.get())
            cisat_db = float(self.CIo_entry.get()) if self.CIo_entry.get().strip() else None
            cisat_lin = 10 ** (cisat_db / 10.0) if cisat_db is not None else None
            bitrate = float(self.bitrate_entry.get()) * 1e6
            demod_loss = float(self.demod_loss_entry.get())
            overhead = float(self.overhead_entry.get())
            uplink_freq_val = _get_optional_float(self.uplink_freq_entry)
            uplink_freq = uplink_freq_val * u.GHz if uplink_freq_val is not None else None
            eirp_gs = _get_optional_float(self.eirp_gs_entry)
            gt_sat = _get_optional_float(self.gt_sat_entry)
            demod_loss_ul = _get_optional_float(self.demod_loss_ul_entry)
            other_att_ul = _get_optional_float(self.other_att_ul_entry)
            uplink_bitrate_val = (
                _get_optional_float(self.uplink_bitrate_entry) if self.uplink_bitrate_entry else None
            )
            uplink_overhead = (
                _get_optional_float(self.uplink_overhead_entry) if self.uplink_overhead_entry else None
            )
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid numerical input: {e}.")
            return

        gs_name = self.gs_var.get()
        lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]
        loss_kwargs = self._pass_loss_settings()
        alt_gs_km = alt_gs_m / 1000.0

        ts = load.timescale()
        tle1 = self.tle1_entry.get().strip()
        tle2 = self.tle2_entry.get().strip()
        try:
            sat = EarthSatellite(tle1, tle2, "SAT", ts)
        except Exception as e:
            messagebox.showerror("TLE Error", f"Invalid TLE data: {e}.")
            return
        gs = wgs84.latlon(lat_gs, lon_gs, elevation_m=alt_gs_m)
        atm_att = atmospheric_attenuation(
            lat_gs,
            lon_gs,
            freq.to(u.GHz).value,
            p,
            d_gs,
            alt_gs_km,
            include_scintillation=self.scint_var.get(),
        )
        uplink_atm_att = None
        if uplink_freq is not None and eirp_gs is not None and gt_sat is not None:
            uplink_atm_att = atmospheric_attenuation(
                lat_gs,
                lon_gs,
                uplink_freq.to(u.GHz).value,
                p,
                d_gs,
                alt_gs_km,
                include_scintillation=self.scint_var.get(),
            )
        self.atm_label_var.set(
            f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}° El: {atm_att:.3f}"
        )
        updated_results = []
        times_dt = self.df_all["Time (UTC)"].tolist()
        sky_times = ts.utc(
            [t.year for t in times_dt],
            [t.month for t in times_dt],
            [t.day for t in times_dt],
            [t.hour for t in times_dt],
            [t.minute for t in times_dt],
            [t.second for t in times_dt],
        )
        altitudes, slant_ranges, dopplers, off_boresight = prepare_topocentric_data(
            sat, gs, sky_times, freq
        )
        for i, t_utc_dt in enumerate(times_dt):
            t_sky = sky_times[i]
            params = calculate_link_budget_parameters(
                t_sky,
                sat,
                gs,
                freq,
                p,
                d_gs,
                alt_gs_km,
                eirp_sat,
                gt_gs,
                demod_loss,
                bitrate,
                overhead,
                cisat_lin,
                other_att,
                atm_att,
                altitudes[i],
                slant_ranges[i],
                dopplers[i],
                off_boresight[i],
                uplink_freq,
                eirp_gs,
                gt_sat,
                demod_loss_ul,
                other_att_ul,
                uplink_atm_att,
                uplink_bitrate_val * 1e6 if uplink_bitrate_val is not None else None,
                uplink_overhead,
                **loss_kwargs,
            )

            rounding_rules = {"Atmospheric Att (dB)": 3, "UL Atmospheric Att (dB)": 3}
            for key in [
                "Elevation (°)",
                "Slant Range (km)",
                "Path Loss (dB)",
                "Atmospheric Att (dB)",
                "S/C Antenna Gain (dBi)",
                "EIRP (dBW)",
                "S/C Antenna Loss (dB)",
                "Off Boresight Angle (°)",
                "Link Losses (dB)",
                "UL Link Losses (dB)",
                "Rx Power (dBW)",
                "C/(No+Io) (dBHz)",
                "Eb/No (dB)",
                "Doppler Shift (kHz)",
                "UL Path Loss (dB)",
                "UL Atmospheric Att (dB)",
                "UL S/C Antenna Loss (dB)",
                "UL Rx Power (dBW)",
                "UL C/No (dBHz)",
                "UL Eb/No (dB)",
            ]:
                if params.get(key) is not None:
                    decimals = rounding_rules.get(key, 2)
                    params[key] = round(params[key], decimals)
            updated_results.append(params)

        self.df_all = pd.DataFrame(updated_results)
        selection = self.contact_listbox.curselection()
        if not selection and len(self.contact_windows) > 0:
            self.contact_listbox.selection_set(0)
            selection = (0,)

        if selection:
            self.on_contact_select(None)
        else:
            self.clear_plot_and_table()

    # ------------------------------------------------------------------
    # Plot / table rendering
    # ------------------------------------------------------------------

    def clear_plot_and_table(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        for widget in self.table_frame.winfo_children():
            widget.destroy()

    def display_contacts(self, contacts):
        self.contact_listbox.delete(0, tk.END)
        filtered_contacts = []
        for contact in contacts:
            if len(contact) == 3:
                start, end, pass_direction = contact
            else:
                start, end = contact
                pass_direction = "A"
            if not self.df_all.empty:
                mask = (self.df_all["Time (UTC)"] >= start) & (self.df_all["Time (UTC)"] <= end)
                df_segment = self.df_all[mask]
                if (
                    not df_segment.empty
                    and df_segment["Elevation (°)"].max() >= MIN_ELEVATION_DEG
                ):
                    filtered_contacts.append((start, end, pass_direction))
        self.contact_windows = filtered_contacts
        for i, (start, end, pass_direction) in enumerate(self.contact_windows):
            mask = (self.df_all["Time (UTC)"] >= start) & (self.df_all["Time (UTC)"] <= end)
            max_elevation = self.df_all.loc[mask, "Elevation (°)"].max()
            minutes = (end - start).total_seconds() / 60.0
            self.contact_listbox.insert(
                tk.END,
                f"Contact {i+1}: {start.strftime('%H:%M:%S')} - {end.strftime('%H:%M:%S')} ({pass_direction})"
                f"  ·  {minutes:.1f} min  ·  max El {max_elevation:.1f}°",
            )

    def on_contact_select(self, event):
        selection = self.contact_listbox.curselection()
        if not selection:
            self.clear_plot_and_table()
            self.clear_uplink_table()
            self.selected_pass_var.set(NO_PASS_SELECTED)
            return

        idx = selection[0]
        start, end, _ = self.contact_windows[idx]
        self.selected_pass_var.set(self.contact_listbox.get(idx))
        if not self.link_budget_enabled:
            return
        mask = (self.df_all["Time (UTC)"] >= start) & (self.df_all["Time (UTC)"] <= end)
        df_pass = self.df_all[mask].copy()
        self.current_table_df = df_pass

        self.clear_uplink_table()
        self.clear_plot_and_table()

        # Constrained layout leaves room for the twin axis of the middle plot,
        # whose right-hand label otherwise ran into the Doppler plot.
        fig, (ax1, ax2, ax3) = plt.subplots(
            1,
            3,
            figsize=(12, 4),
            gridspec_kw={"width_ratios": [1, 1, 1]},
            layout="constrained",
        )

        tk_bg_color = self.root.cget("bg")
        try:
            r, g, b = self.root.winfo_rgb(tk_bg_color)
            mpl_bg_color = f"#{r // 256:02x}{g // 256:02x}{b // 256:02x}"
        except tk.TclError:
            mpl_bg_color = "#F0F0F0"
        fig.patch.set_facecolor(mpl_bg_color)

        ax1.plot(df_pass["Time (UTC)"], df_pass["Elevation (°)"], label="Elevation (°)", color=MATLAB_COLORS["blue"])
        ax1.set_ylabel("Elevation (°)", color=MATLAB_COLORS["blue"], fontsize=8)
        ax1.tick_params(axis="both", labelcolor=MATLAB_COLORS["blue"], labelsize=8)
        ax1.set_xlabel("Time (UTC)", fontsize=8)
        ax1.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.grid(True)
        ax1.set_title("Elevation Angle")

        ax2_1 = ax2
        ax2_2 = ax2.twinx()
        ax2_1.plot(df_pass["Time (UTC)"], df_pass["Eb/No (dB)"], label="Eb/No (dB)", color=MATLAB_COLORS["orange"])
        ax2_1.set_ylabel("Eb/No (dB)", color=MATLAB_COLORS["orange"], fontsize=8)
        ax2_1.tick_params(axis="both", labelcolor=MATLAB_COLORS["orange"], labelsize=8)
        ax2_2.plot(
            df_pass["Time (UTC)"],
            df_pass["C/(No+Io) (dBHz)"],
            label="C/(No+Io) (dBHz)",
            color=MATLAB_COLORS["green"],
            linestyle="--",
        )
        ax2_2.set_ylabel("C/(No+Io) (dBHz)", color=MATLAB_COLORS["green"], fontsize=8)
        ax2_2.tick_params(axis="both", labelcolor=MATLAB_COLORS["green"], labelsize=8)
        ax2.set_xlabel("Time (UTC)", fontsize=8)
        ax2.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.grid(True)
        ax2.set_title("Link Quality (Eb/No & C/(No+Io))")
        lines, labels = ax2_1.get_legend_handles_labels()
        lines2, labels2 = ax2_2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc="best")

        ax3.plot(df_pass["Time (UTC)"], df_pass["Doppler Shift (kHz)"], color=MATLAB_COLORS["purple"])
        ax3.set_xlabel("Time (UTC)", fontsize=8)
        ax3.set_ylabel("Doppler (kHz)", color=MATLAB_COLORS["purple"], fontsize=8)
        ax3.tick_params(axis="both", labelcolor=MATLAB_COLORS["purple"], labelsize=8)
        ax3.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax3.grid(True)
        ax3.set_title("Doppler Shift")
        # Rotated time labels, set on the axes: fig.autofmt_xdate() would call
        # subplots_adjust, which the constrained layout does not allow.
        for ax in (ax1, ax2, ax3):
            ax.tick_params(axis="x", labelrotation=30)

        canvas_plot_widget = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas_plot_widget.draw()
        canvas_plot_widget.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)

        display_columns = [
            "Time (UTC)",
            "Elevation (°)",
            "Slant Range (km)",
            "Path Loss (dB)",
            "Atmospheric Att (dB)",
            "Off Boresight Angle (°)",
            "S/C Antenna Gain (dBi)",
            "EIRP (dBW)",
            "Link Losses (dB)",
            "Rx Power (dBW)",
            "C/(No+Io) (dBHz)",
            "Eb/No (dB)",
            "Doppler Shift (kHz)",
        ]
        table = ttk.Treeview(self.table_frame, columns=display_columns, show="headings")
        table.tag_configure("evenrow", background=PALETTE["row_even"])
        table.tag_configure("oddrow", background=PALETTE["row_odd"])
        for col in display_columns:
            table.heading(col, text=col)
            if col == "Time (UTC)":
                table.column(col, anchor="center", width=120)
            else:
                table.column(col, anchor="center", width=100)
        for i, (_, row) in enumerate(df_pass.iterrows()):
            formatted_time = row["Time (UTC)"].strftime("%Y-%m-%d %H:%M:%S")
            values = [formatted_time] + [row.get(col, "") for col in display_columns[1:]]
            table.insert("", "end", values=values, tags=("evenrow" if i % 2 == 0 else "oddrow",))
        table_vscroll = ttk.Scrollbar(self.table_frame, orient=tk.VERTICAL, command=table.yview)
        table_hscroll = ttk.Scrollbar(self.table_frame, orient=tk.HORIZONTAL, command=table.xview)
        table.configure(yscrollcommand=table_vscroll.set, xscrollcommand=table_hscroll.set)
        table_vscroll.pack(side=tk.RIGHT, fill=tk.Y)
        table_hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        table.pack(fill=tk.BOTH, expand=True)

    # ------------------------------------------------------------------
    # Miscellaneous actions
    # ------------------------------------------------------------------

    def show_about(self):
        """Show the version, the author and the licence of the program."""
        lines = [
            f"Satellite Link Budget Tool {APP_VERSION}",
            "Orbit prediction & downlink/uplink budget analysis",
            "",
            f"Author: {APP_AUTHOR}",
        ]
        license_text = self.license_status.summary()
        if license_text:
            lines += ["", license_text]
        messagebox.showinfo("About", "\n".join(lines), parent=self.root)

    def exit_app(self):
        """Close the program: the Exit button, Ctrl+Q and the window close button."""

        plt.close("all")
        self.root.quit()
        self.root.destroy()
        sys.exit()

    def show_antenna_pattern(self):
        popup = tk.Toplevel(self.root)
        popup.title("Satellite Antenna Pattern")
        popup.geometry("600x400")
        popup.configure(bg=PALETTE["bg"])

        from calculations import ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS, marker="o", linestyle="-")
        ax.set_xlabel("Angle from nadir (boresight) [deg]", fontsize=8)
        ax.set_ylabel("Gain [dBi]", fontsize=8)
        ax.set_title("Satellite Antenna Gain Pattern")
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True)
        canvas_popup = FigureCanvasTkAgg(fig, master=popup)
        canvas_popup.draw()
        canvas_popup.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)

    def load_antenna_pattern_file(self):
        """Load the satellite antenna pattern (angle from nadir, gain in dBi) from a CSV file."""
        file_path = filedialog.askopenfilename(
            title="Select Satellite Antenna Pattern File",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")],
        )
        if not file_path:
            return
        try:
            load_antenna_pattern(file_path)
            self.pattern_file_var.set(f"Pattern: {os.path.abspath(file_path)}")
            messagebox.showinfo("Success", "Satellite antenna pattern loaded successfully.")
            self.set_analysis_stale()
            # The static tabs read the Tx antenna gain from the pattern.
            for refresh in self._static_lb_refreshers:
                refresh()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load antenna pattern: {e}")

    def export_table_csv(self):
        """Export the currently displayed table to a CSV file."""
        if self.current_table_df.empty:
            messagebox.showwarning("Warning", "Please select a contact window first.")
            return
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                project = self.project_var.get().strip()
                if project:
                    # Project and pass on the first lines, then the table.
                    f.write(f"Project,{project}\n")
                    f.write(f"Pass,{self.selected_pass_var.get()}\n\n")
                self.current_table_df.to_csv(f, index=False)
            messagebox.showinfo("Export Successful", f"Data exported to {file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export CSV: {e}")

    # ------------------------------------------------------------------
    # Fixed-elevation preliminary link budget (no TLE required)
    # ------------------------------------------------------------------

    def show_fixed_downlink_link_budget(self):
        """Show the static (fixed-elevation) downlink link budget tab."""

        self.main_tabs.select(self.static_dl_tab)

    def show_fixed_uplink_link_budget(self):
        """Show the static (fixed-elevation) uplink link budget tab."""

        self.main_tabs.select(self.static_ul_tab)

    def _build_static_downlink_view(self, parent):
        """The fixed-elevation preliminary downlink link budget, in ``parent``."""

        self._build_fixed_link_budget_view(
            parent,
            "Fixed-Elevation Downlink Link Budget",
            freq_entry=self.freq_entry,
            freq_label="Frequency [GHz]",
            eirp_entry=self.eirp_sat_entry,
            eirp_label="EIRP SAT [dBW]",
            gt_entry=self.gt_gs_entry,
            gt_label="G/T GS [dB/K]",
            demod_loss_entry=self.demod_loss_entry,
            demod_loss_label="Demodulator Loss [dB]",
            bitrate_entry=self.bitrate_entry,
            bitrate_label="Coded Bit Rate [Mbps]",
            overhead_entry=self.overhead_entry,
            overhead_label="Overhead (coded / info)",
            other_att_entry=self.other_att_entry,
            other_att_label="Other Attenuations [dB]",
            rolloff_entry=self.rolloff_entry,
            spectral_eff_entry=self.spectral_eff_entry,
            rx_is_ground_station=True,
        )

    def _build_static_uplink_view(self, parent):
        """The fixed-elevation preliminary uplink link budget, in ``parent``."""

        self._build_fixed_link_budget_view(
            parent,
            "Fixed-Elevation Uplink Link Budget",
            freq_entry=self.uplink_freq_entry,
            freq_label="Uplink Frequency [GHz]",
            eirp_entry=self.eirp_gs_entry,
            eirp_label="EIRP GS [dBW]",
            gt_entry=self.gt_sat_entry,
            gt_label="G/T SAT [dB/K]",
            demod_loss_entry=self.demod_loss_ul_entry,
            demod_loss_label="Sat Demod Loss [dB]",
            bitrate_entry=self.uplink_bitrate_entry,
            bitrate_label="UL Coded Bit Rate [Mbps]",
            overhead_entry=self.uplink_overhead_entry,
            overhead_label="UL Overhead (coded / info)",
            other_att_entry=self.other_att_ul_entry,
            other_att_label="Uplink Other Attenuations [dB]",
            rolloff_entry=self.uplink_rolloff_entry,
            spectral_eff_entry=self.uplink_spectral_eff_entry,
            rx_is_ground_station=False,
        )

    def _build_fixed_link_budget_view(
        self,
        parent: tk.Widget,
        title: str,
        freq_entry: ttk.Entry,
        freq_label: str,
        eirp_entry: ttk.Entry,
        eirp_label: str,
        gt_entry: ttk.Entry,
        gt_label: str,
        demod_loss_entry: ttk.Entry,
        demod_loss_label: str,
        bitrate_entry: ttk.Entry,
        bitrate_label: str,
        overhead_entry: ttk.Entry,
        overhead_label: str,
        other_att_entry: ttk.Entry,
        other_att_label: str,
        rolloff_entry: ttk.Entry,
        spectral_eff_entry: ttk.Entry,
        rx_is_ground_station: bool = True,
    ):
        """Build, in ``parent`` (a main window tab), a TLE-independent fixed-elevation link budget.

        This is the standard preliminary "worst case" check: pick a single
        (typically minimum) elevation angle and satellite altitude, derive
        the slant range purely from that geometry, and verify a link margin
        exists before running the full TLE-based pass analysis. Frequency,
        G/T, bit rate and the other panel values are shared with the pass
        analysis tab of the same direction so they don't need to be retyped.

        The view's own inputs are grouped as in a formal static link budget
        sheet (geometry, data channel and requirements, transmitter,
        propagation and PFD limit, receiver) plus favourable/adverse
        tolerances for the statistical margins. Each section, and each part
        of the results (summary, link budget, PFD check, link budget with
        tolerances), is a page picked from a sidebar; the action buttons and
        the status stay visible above the pages. The inputs can be saved
        into / loaded from a parameters JSON file (as the
        ``fixed_elevation_downlink`` / ``fixed_elevation_uplink`` section),
        and are filled from ``parameters.json`` at startup and whenever a
        parameters file is imported.
        """

        section = calculations.FIXED_LINK_BUDGET_SECTIONS["downlink" if rx_is_ground_station else "uplink"]
        tx_side, rx_side = ("S/C", "GS") if rx_is_ground_station else ("GS", "S/C")

        # The view lives in a tab of the main window for the whole session.
        win = parent

        recap_font_family = _pick_font(["Segoe UI", "Helvetica Neue", "Helvetica", "Arial"])
        direction_text = "Downlink (S/C → GS)" if rx_is_ground_station else "Uplink (GS → S/C)"

        # Layout: an action bar (buttons and status), then a sidebar listing
        # the input and result pages next to the page shown.
        action_bar = ttk.Frame(win, style="Toolbar.TFrame", padding=(16, 10))
        action_bar.pack(fill=tk.X)
        button_row = ttk.Frame(action_bar, style="Toolbar.TFrame")
        button_row.pack(side=tk.LEFT)
        status_box = ttk.Frame(action_bar, style="Toolbar.TFrame")
        status_box.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(16, 0))
        tk.Frame(win, bg=PALETTE["border"], height=1).pack(fill=tk.X)

        main_area = ttk.Frame(win)
        main_area.pack(fill=tk.BOTH, expand=True)
        sidebar = tk.Frame(main_area, bg=PALETTE["surface"], width=236)
        sidebar.pack(side=tk.LEFT, fill=tk.Y)
        sidebar.pack_propagate(False)
        tk.Frame(main_area, bg=PALETTE["border"], width=1).pack(side=tk.LEFT, fill=tk.Y)
        pages_book = ttk.Notebook(main_area, style="Pages.TNotebook")
        pages_book.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # key -> {"container", "body", "canvas", "nav", "label", "badge", "enabled"}
        pages: dict[str, dict] = {}
        current_page = {"key": None}
        nav_font = (recap_font_family, 10)
        nav_font_selected = (recap_font_family, 10, "bold")

        def _paint_nav(key: str, hover: bool = False) -> None:
            page = pages[key]
            if not page["enabled"]:
                bg, fg, bar, font = PALETTE["surface"], PALETTE["nav_disabled"], PALETTE["surface"], nav_font
            elif current_page["key"] == key:
                bg, fg, bar, font = PALETTE["accent_light"], PALETTE["accent_dark"], PALETTE["accent"], nav_font_selected
            else:
                bg = PALETTE["surface_alt"] if hover else PALETTE["surface"]
                fg, bar, font = PALETTE["text"], bg, nav_font
            page["nav"].configure(bg=bg)
            page["bar"].configure(bg=bar)
            page["label"].configure(bg=bg, fg=fg, font=font)
            page["badge"].configure(bg=bg)

        def _show_page(key: str) -> None:
            if key not in pages or not pages[key]["enabled"]:
                return
            previous = current_page["key"]
            current_page["key"] = key
            pages_book.select(pages[key]["container"])
            pages[key]["canvas"].yview_moveto(0)
            for other in {previous, key} - {None}:
                _paint_nav(other)

        def _set_page_enabled(key: str, enabled: bool) -> None:
            pages[key]["enabled"] = enabled
            pages[key]["nav"].configure(cursor="hand2" if enabled else "arrow")
            _paint_nav(key)

        def _set_page_badge(key: str, text: str = "", colour: str | None = None) -> None:
            pages[key]["badge"].configure(text=text, fg=colour or PALETTE["text_muted"])

        def _nav_caption(text: str) -> None:
            tk.Label(
                sidebar,
                text=text,
                bg=PALETTE["surface"],
                fg=PALETTE["text_muted"],
                font=(recap_font_family, 8, "bold"),
                anchor="w",
            ).pack(fill=tk.X, padx=18, pady=(16, 4))

        def _add_page(key: str, nav_text: str, heading: str = "", description: str = "") -> ttk.Frame:
            """Add a page and its sidebar entry; return the frame for its content."""

            container, body, canvas = self._make_scrollable_frame(pages_book)
            pages_book.add(container)
            nav = tk.Frame(sidebar, bg=PALETTE["surface"], cursor="hand2")
            nav.pack(fill=tk.X)
            bar = tk.Frame(nav, bg=PALETTE["surface"], width=4)
            bar.pack(side=tk.LEFT, fill=tk.Y)
            label = tk.Label(nav, text=nav_text, bg=PALETTE["surface"], anchor="w", padx=14, pady=7)
            label.pack(side=tk.LEFT, fill=tk.X, expand=True)
            badge = tk.Label(nav, text="", bg=PALETTE["surface"], font=(recap_font_family, 10, "bold"), padx=12)
            badge.pack(side=tk.RIGHT)
            pages[key] = {
                "container": container, "body": body, "canvas": canvas, "nav": nav, "bar": bar,
                "label": label, "badge": badge, "enabled": True,
            }
            for widget in (nav, bar, label, badge):
                widget.bind("<Button-1>", lambda _e, k=key: _show_page(k))
                widget.bind("<Enter>", lambda _e, k=key: k != current_page["key"] and _paint_nav(k, hover=True))
                widget.bind("<Leave>", lambda _e, k=key: _paint_nav(k))
            _paint_nav(key)

            content = ttk.Frame(body, padding=(22, 18, 22, 22))
            content.pack(fill=tk.BOTH, expand=True)
            pages[key]["content"] = content
            if heading:
                ttk.Label(content, text=heading, style="PageTitle.TLabel").pack(anchor="w")
            if description:
                ttk.Label(
                    content, text=description, style="PageSub.TLabel", wraplength=760, justify="left"
                ).pack(anchor="w", pady=(2, 0))
            if heading or description:
                ttk.Frame(content, height=12).pack(fill=tk.X)
            return content

        def _page_of(widget) -> str | None:
            """The key of the page that contains ``widget``, if any."""

            containers = {str(page["container"]): key for key, page in pages.items()}
            while widget is not None:
                if str(widget) in containers:
                    return containers[str(widget)]
                widget = widget.master
            return None

        # Mouse-wheel scrolling goes to the page shown while this tab is selected.
        self._wheel_targets[str(win)] = (
            lambda: pages[current_page["key"]]["canvas"] if current_page["key"] else None
        )

        # Sections in the order of a static link budget sheet. Each item is
        # one of:
        #   ("field", key, label, default text, placeholder, tooltip) -- an
        #       input of this view. Fields with a real default ("0", the
        #       minimum elevation, ...) start filled in; fields whose empty
        #       value means "not used" (or "derive it automatically") show a
        #       greyed-out example that is never treated as real input.
        #   ("panel", label, panel entry, tooltip) -- a parameters panel
        #       field, edited through the panel's own StringVar, with the
        #       panel's handling (derived values, stale flags) on every edit.
        #   ("station",) / ("scintillation",) -- the panel's ground station
        #       selector and scintillation switch.
        def F(key, label, default, placeholder, tooltip):
            return ("field", key, label, default, placeholder, tooltip)

        def L(label, panel_entry, tooltip=""):
            return ("panel", label, panel_entry, tooltip)

        field_groups = [
            (
                "Ground Station & Geometry",
                [
                    ("station",),
                    F(
                        "elevation_deg", "Elevation Angle [deg]", f"{MIN_ELEVATION_DEG:g}", None,
                        "Required. Elevation of the satellite seen from the ground station; the "
                        "minimum operational elevation gives the worst case (longest slant range).",
                    ),
                    F(
                        "sat_altitude_km", "Satellite Altitude [km]", "", "e.g. 628",
                        "Required. Circular-orbit altitude used, with the elevation, to derive the "
                        "slant range.",
                    ),
                    L(
                        "Antenna Diameter GS [m]", self.d_gs_entry,
                        "Used by the ITU-R atmospheric model"
                        + (" and for the default Rx beamwidth." if rx_is_ground_station else "."),
                    ),
                    L(
                        "Link Availability [%]", self.LA_entry,
                        "Weather availability: rain/cloud attenuation is computed for "
                        "p = 100 − availability [% of time].",
                    ),
                    ("scintillation",),
                ],
            ),
            (
                "Data Channel & Requirements",
                (
                    # CCSDS TM framing applies to the downlink (telemetry) only.
                    [
                        ("framing_mode",),
                        F(
                            "ccsds_interleaving", "RS Interleaving Depth I", "5", None,
                            "RS(255,223) interleaving depth from the mission ICD: the transfer frame is "
                            "223·I bytes (I = 5 → 1115 bytes). Used in CCSDS mode.",
                        ),
                        ("ccsds_options",),
                        F(
                            "ccsds_packet_length", "Source Packet Length [B]", "", "empty = not deducted",
                            "Total length of the Source Packets (header to error control); packets span "
                            "transfer frames. Each packet's 6-byte header, data field header and packet "
                            "error control are not user data: they lower the information bit rate, not "
                            "the margin (Eb/No is referred to the frame bit rate). Leave empty not to "
                            "deduct packet overheads.",
                        ),
                        F(
                            "ccsds_packet_sec_header", "Data Field Header (CUC) [B]", "0", None,
                            "Length of the packet data field header (e.g. a CUC time code: optional "
                            "1-byte P-field + coarse and fine time octets) when the Data Field Header "
                            "Flag is 1. Used with the packet length.",
                        ),
                        ("framing_used",),
                    ]
                    if rx_is_ground_station
                    else []
                )
                + [
                    L(
                        bitrate_label, bitrate_entry,
                        "Coded (channel) bit rate of the transmitter, channel coding included. In "
                        "CCSDS mode the transfer frame and information bit rates are derived from it.",
                    ),
                    L(
                        overhead_label, overhead_entry,
                        "Coded bit rate / information (user data) bit rate = formatting factor × "
                        "channel coding factor. Computed in CCSDS mode.",
                    ),
                    F(
                        "formatting_overhead", "Formatting Factor (frame / user data)", "", "e.g. 1.02",
                        "Transfer frame bit rate divided by the information (user data) bit rate: "
                        "frame header and trailer, and packet overheads. Eb/No is referred to the "
                        "transfer frame bit rate = information bit rate × this factor. Empty = 1. "
                        "Computed in CCSDS mode.",
                    ),
                    ("computed", "coding_factor", "Coding Factor (coded / frame)"),
                    ("computed", "frame_rate", "Transfer Frame Bit Rate [kbit/s]"),
                    ("computed", "info_rate", "Information Bit Rate [kbit/s]"),
                ]
                # Roll-off and spectral efficiency only set the occupied bandwidth
                # for the PFD check, which applies to the downlink only.
                + (
                    [
                        L(
                            "Roll-off Factor", rolloff_entry,
                            "With the spectral efficiency, gives the occupied bandwidth used "
                            "for the PFD in 4 kHz.",
                        ),
                        L(
                            "Spectral Efficiency [bps/Hz]", spectral_eff_entry,
                            "Coded bits per second per Hz of the modulation and filtering; with the "
                            "roll-off gives the occupied bandwidth used for the PFD in 4 kHz.",
                        ),
                    ]
                    if rx_is_ground_station
                    else []
                )
                + [
                    F(
                        "required_ebno_db", "Required Eb/No [dB]", "", "e.g. 2.6",
                        "Eb/No threshold of the modulation/coding for the required error rate. "
                        "When set, the Data Recovery Margin (received minus required Eb/No) is "
                        "computed.",
                    ),
                    F(
                        "required_margin_db", "Required Nominal Margin [dB]", "3", None,
                        "Requirement on the nominal Data Recovery Margin (typically ≥ 3 dB): the "
                        "margin is shown green/red against it. MEAN − 3σ and Worst Case RSS are "
                        "checked against 0 dB.",
                    ),
                ],
            ),
            (
                f"Transmitter ({tx_side})",
                [
                    ("eirp_mode",),
                    L(
                        freq_label, freq_entry,
                        "Carrier frequency: sets the free space loss and the atmospheric attenuation.",
                    ),
                    L(
                        eirp_label, eirp_entry,
                        "Used when the EIRP source is 'EIRP value'.",
                    ),
                    F(
                        "tx_power_w", "Transmitter Power [W]", "", "e.g. 3.3",
                        "Tx chain (EIRP source 'Tx chain'): EIRP = power − circuit loss − VSWR "
                        "loss + antenna gain.",
                    ),
                    F(
                        "ant_circuit_loss_db", "Antenna Circuit Loss [dB]", "0", None,
                        "Losses between transmitter and antenna (cables, filters, couplers). "
                        "Part of the Tx chain.",
                    ),
                    F(
                        "vswr", "VSWR (:1)", "", "e.g. 1.5",
                        "Overall VSWR, e.g. 1.5 for 1.5:1; gives the mismatch loss "
                        "−10·log10(1 − Γ²). Part of the Tx chain.",
                    ),
                    F(
                        "ant_gain_dbi", "Antenna Gain [dBi]", "", "e.g. 0",
                        "Tx antenna gain towards the receiver. Part of the Tx chain (required "
                        "there)."
                        + (
                            " Can be read from the satellite antenna pattern instead (option below)."
                            if rx_is_ground_station
                            else ""
                        ),
                    ),
                ]
                # The satellite antenna pattern gives the Tx gain on the downlink only.
                + ([("gain_source",)] if rx_is_ground_station else [])
                + [
                    ("eirp_used",),
                    F(
                        "ant_axial_ratio_db", "Tx Antenna Axial Ratio [dB]", "", "e.g. 3",
                        "With the Rx axial ratio, gives the polarisation mismatch loss when that "
                        "field is left empty.",
                    ),
                    F(
                        "tx_pointing_loss_db", "Tx Antenna Pointing Loss [dB]", "0", None,
                        "Transmit antenna pointing / mispointing loss, subtracted from Rx Power "
                        "(0 for an omnidirectional antenna).",
                    ),
                ],
            ),
            (
                "Propagation & PFD Limit" if rx_is_ground_station else "Propagation",
                [
                    ("iono_mode",),
                    F(
                        "ionospheric_loss_db", "Ionospheric Loss [dB]", "0", None,
                        "Manual ionospheric allowance (scintillation and absorption), subtracted "
                        "from Rx Power. Used when the ionospheric loss is 'Manual'.",
                    ),
                    F(
                        "iono_s4_ref", "Ionospheric S4 (reference)", "0.3", None,
                        "S4 scintillation index exceeded for the wanted percentage of time, at the "
                        "reference frequency (published statistics are mostly at GPS L1). 0.3 is the "
                        "99.9% value measured in Kiruna (auroral oval, active solar years) used by "
                        "ESA AWS-TN-ESA-SYS-0052; it is a conservative upper bound for Svalbard. "
                        "Used when the ionospheric loss is 'From S4'.",
                    ),
                    F(
                        "iono_s4_ref_freq_mhz", "S4 Reference Frequency [MHz]", f"{calculations.GPS_L1_MHZ:g}", None,
                        "Frequency at which the reference S4 was measured (GPS L1 = 1575.42 MHz). "
                        "S4 is scaled to the link frequency as f^-1.5 (ITU-R P.531).",
                    ),
                    ("iono_used",),
                    L(
                        other_att_label, other_att_entry,
                        "Any other static attenuation, subtracted from Rx Power. Do not put here a "
                        "total atmospheric allowance (e.g. 3 dB) if the atmospheric, ionospheric "
                        "and polarisation losses are also computed: they would be counted twice.",
                    ),
                    F(
                        "polarisation_loss_db", "Polarisation Mismatch Loss [dB]", "", "auto from axial ratios",
                        "Subtracted from Rx Power. Leave empty to compute it from the Tx and Rx "
                        "axial ratios (worst-case orientation of the polarisation ellipses, with the "
                        "antenna gains taken as co-polar/circular gains, as in the reference "
                        "sheet), or 0 when those are not both set.",
                    ),
                    F(
                        "multipath_loss_db", "Multipath Losses [dB]", "0", None,
                        "Allowance for multipath/ground reflections (terrain, snow, nearby "
                        "structures), mainly at low elevation. Subtracted from Rx Power. It is "
                        "site-specific: simulate it for the actual antenna site, with and without "
                        "snow, and confirm it by measurement (ESA AWS-TN-ESA-SYS-0052).",
                    ),
                ]
                + (
                    [
                        F(
                            "pfd_limit", "PFD Limit [dBW/m²/4kHz]", "", "e.g. -154",
                            "Regulatory power flux density limit in a 4 kHz reference bandwidth "
                            "(ITU Radio Regulations Art. 21, ECSS-E-ST-50-05C). It depends on the "
                            "frequency band and on the elevation angle. When set, the PFD Margin "
                            "is computed.",
                        ),
                    ]
                    # PFD limits apply to space stations (downlink) only.
                    if rx_is_ground_station
                    else []
                ),
            ),
            (
                f"Receiver ({rx_side})",
                [
                    L(
                        gt_label, gt_entry,
                        "Receiving figure of merit: antenna gain over system noise temperature."
                        + (
                            ""
                            if rx_is_ground_station
                            else " With the option below it is the value at the satellite antenna "
                            "boresight (pattern peak)."
                        ),
                    ),
                ]
                # The satellite antenna pattern gives the Rx gain on the uplink only.
                + ([] if rx_is_ground_station else [("gain_source",)])
                + [
                    L(
                        demod_loss_label, demod_loss_entry,
                        "Receiver (demodulator) implementation loss. Lowers Eb/No.",
                    ),
                    F(
                        "modulation_degradation_db", "Modulation Degradation [dB]", "0", None,
                        "Modulator/filter implementation loss, added to the demodulator loss. "
                        "Lowers Eb/No.",
                    ),
                    F(
                        "rx_axial_ratio_db", "Rx Antenna Axial Ratio [dB]", "", "e.g. 1",
                        "With the Tx axial ratio, gives the polarisation mismatch loss when that "
                        "field is left empty.",
                    ),
                    F(
                        "rx_depointing_deg", "Rx Antenna Depointing Angle [deg]", "", "e.g. 0.1",
                        "Tracking/pointing error of the receive antenna; gives the pointing loss "
                        "12·(angle / beamwidth)², subtracted from Rx Power.",
                    ),
                    F(
                        "rx_beamwidth_deg", "Rx Antenna 3dB Beamwidth [deg]", "",
                        "auto from GS diameter" if rx_is_ground_station else "e.g. 10",
                        "Used together with the depointing angle. "
                        + (
                            "Leave empty to use 70·λ/D from the GS antenna diameter."
                            if rx_is_ground_station
                            else "Required when the depointing angle is set."
                        ),
                    ),
                ],
            ),
        ]

        # (page key, sidebar name, description) of each section, in the order of field_groups.
        input_pages = [
            (
                "geometry",
                "Station & Geometry",
                "Where the link is evaluated: the ground station, the fixed elevation and the satellite "
                "altitude (together they give the slant range), and the weather availability used for "
                "rain and cloud attenuation.",
            ),
            (
                "data",
                "Data Channel",
                "Bit rate, framing and overheads of the data channel, and the Eb/No and margin the link "
                "has to achieve.",
            ),
            (
                "transmitter",
                f"Transmitter ({tx_side})",
                "EIRP, entered directly or computed from the transmitter chain (power, losses, antenna gain).",
            ),
            (
                "propagation",
                "Propagation & PFD" if rx_is_ground_station else "Propagation",
                "Losses between the two antennas besides free space: atmosphere (ITU-R), ionosphere, "
                "polarisation, multipath and other attenuations"
                + (", and the regulatory PFD limit." if rx_is_ground_station else "."),
            ),
            (
                "receiver",
                f"Receiver ({rx_side})",
                "Figure of merit of the receiving side and the losses of the receiver and demodulator.",
            ),
        ]
        entries: dict[str, ttk.Entry] = {}
        entry_labels: dict[str, str] = {}
        entry_defaults: dict[str, str] = {}
        panel_specs: list[tuple[str, ttk.Entry]] = []
        # key -> (label widget, entry widget, label colour), to grey out the
        # fields the calculation currently ignores. Panel fields use the key
        # "panel:<label>".
        field_widgets: dict[str, tuple[ttk.Label, ttk.Entry, str]] = {}
        # key -> (info icon, its default text), for fields whose icon can show
        # a live breakdown of the computed value.
        info_icons: dict[str, tuple[ttk.Label, str]] = {}
        # Read-only boxes with values derived from the bit rate and overheads.
        computed_boxes: dict[str, tuple[tk.StringVar, ttk.Label]] = {}
        eirp_mode_var = tk.StringVar(value="eirp")
        gain_pattern_var = tk.BooleanVar(value=False)
        gain_pattern_text_var = tk.StringVar(value="")
        # The station of this view is its own (saved with its inputs), so
        # changing it here never touches the pass analysis.
        station_var = tk.StringVar(value=self.gs_var.get())
        iono_mode_var = tk.StringVar(value="manual")
        framing_mode_var = tk.StringVar(value="manual")
        ccsds_ocf_var = tk.BooleanVar(value=True)
        ccsds_fecf_var = tk.BooleanVar(value=False)
        ccsds_pec_var = tk.BooleanVar(value=True)
        framing_used_var = tk.StringVar(value="")
        ccsds_option_buttons: list = []
        iono_rows: dict[str, tuple[list, tk.StringVar]] = {}
        # Tropospheric scintillation is also this view's own choice: the
        # ITU-R P.618 model does not apply to every link (see the results notes).
        scint_var = tk.BooleanVar(value=self.scint_var.get())
        eirp_used_var = tk.StringVar(value="")
        passes_tab = "Downlink" if rx_is_ground_station else "Uplink"
        shared_note = "Shared with the pass analysis tabs."
        _nav_caption("INPUTS")
        for group_index, (group_title, items) in enumerate(field_groups):
            page_key, nav_text, page_description = input_pages[group_index]
            page_content = _add_page(page_key, nav_text, group_title, page_description)
            group_frame = ttk.Frame(page_content)
            group_frame.pack(fill=tk.X)
            group_frame.grid_columnconfigure(1, weight=1)
            group_frame.grid_columnconfigure(3, weight=1)
            r, c = 0, 0
            for item in items:
                kind = item[0]
                if kind in (
                    "eirp_mode", "eirp_used", "iono_mode", "iono_used", "framing_mode", "ccsds_options",
                    "framing_used", "gain_source",
                ):
                    # Full-width rows.
                    if c:
                        r, c = r + 1, 0
                    if kind == "framing_mode":
                        framing_row = ttk.Frame(group_frame)
                        framing_row.grid(row=r, column=0, columnspan=4, sticky="w", pady=(0, 6))
                        ttk.Label(
                            framing_row, text="Data framing:", font=(recap_font_family, 9, "bold")
                        ).pack(side=tk.LEFT, padx=(0, 8))
                        for value, text in (
                            ("manual", "Manual (bit rate and overheads below)"),
                            ("ccsds", "CCSDS TM: RS(255,223) + convolutional 1/2"),
                        ):
                            ttk.Radiobutton(
                                framing_row,
                                text=text,
                                value=value,
                                variable=framing_mode_var,
                                command=lambda: _refresh_states(),
                            ).pack(side=tk.LEFT, padx=(0, 12))
                        _info_icon(
                            framing_row,
                            "CCSDS TM: from the coded (transmitter) bit rate and the RS interleaving "
                            "depth I, the transfer frame structure (223·I bytes: 6-byte primary "
                            "header, optional OCF 4 bytes and FECF 2 bytes, data field) and the "
                            "channel coding (32-bit ASM per RS(255,223) codeblock of 255·I bytes, "
                            "convolutional rate 1/2) give the formatting overhead and the total "
                            f"overhead, filled in below (the overhead also in the {passes_tab} tab), and "
                            "the transfer frame (Eb/No reference) and information bit rates.",
                        ).pack(side=tk.LEFT)
                    elif kind == "ccsds_options":
                        options_row = ttk.Frame(group_frame)
                        options_row.grid(row=r, column=0, columnspan=4, sticky="w", pady=2)
                        for text, var in (
                            ("Frame OCF (4 B)", ccsds_ocf_var),
                            ("Frame FECF (2 B)", ccsds_fecf_var),
                            ("Packet Error Control (2 B)", ccsds_pec_var),
                        ):
                            button = ttk.Checkbutton(
                                options_row, text=text, variable=var, command=lambda: _refresh_states()
                            )
                            button.pack(side=tk.LEFT, padx=(0, 16))
                            ccsds_option_buttons.append(button)
                        _info_icon(
                            options_row,
                            "As in the mission ICD. Frame OCF and FECF: the optional transfer frame "
                            "trailer (the OCF carries the CLCW; the FECF is often omitted with RS "
                            "coding). Packet Error Control: the 16-bit CRC at the end of each Source "
                            "Packet, deducted when the packet length is given.",
                        ).pack(side=tk.LEFT)
                    elif kind == "framing_used":
                        ttk.Label(
                            group_frame,
                            textvariable=framing_used_var,
                            font=(recap_font_family, 9, "bold"),
                            foreground=PALETTE["accent_dark"],
                            wraplength=760,
                            justify="left",
                        ).grid(row=r, column=0, columnspan=4, sticky="w", pady=(4, 4))
                    elif kind == "iono_mode":
                        iono_row = ttk.Frame(group_frame)
                        iono_row.grid(row=r, column=0, columnspan=4, sticky="w", pady=(0, 6))
                        ttk.Label(
                            iono_row, text="Ionospheric loss:", font=(recap_font_family, 9, "bold")
                        ).pack(side=tk.LEFT, padx=(0, 8))
                        for value, text in (
                            ("manual", "Manual"),
                            ("s4", "From S4 scintillation index (ITU-R P.531)"),
                        ):
                            ttk.Radiobutton(
                                iono_row,
                                text=text,
                                value=value,
                                variable=iono_mode_var,
                                command=lambda: _refresh_states(),
                            ).pack(side=tk.LEFT, padx=(0, 12))
                        _info_icon(
                            iono_row,
                            "Manual: the Ionospheric Loss value below. From S4: the ionospheric "
                            "scintillation fade derived from an S4 index as in ESA "
                            "AWS-TN-ESA-SYS-0052: S4 scaled to the link frequency as f^-1.5, "
                            "peak-to-peak fluctuation 27.5·S4^1.26 dB, fade = half of it. The "
                            "station's geomagnetic zone is shown below as guidance for S4.",
                        ).pack(side=tk.LEFT)
                    elif kind == "iono_used":
                        # Small summary table: one row per piece of information.
                        iono_box = ttk.Frame(group_frame)
                        iono_box.grid(row=r, column=0, columnspan=4, sticky="w", pady=(4, 6))
                        for i, (row_key, row_label) in enumerate(
                            (
                                ("zone", "Station zone"),
                                ("guidance", "S4 guidance"),
                                ("s4", "S4 at link frequency"),
                                ("fluctuation", "Peak-to-peak fluctuation"),
                                ("loss", "Ionospheric loss used"),
                            )
                        ):
                            emphasis = row_key == "loss"
                            name_label = ttk.Label(
                                iono_box,
                                text=f"{row_label}:",
                                foreground=PALETTE["text"] if emphasis else PALETTE["text_muted"],
                                font=(recap_font_family, 9, "bold") if emphasis else (recap_font_family, 9),
                            )
                            name_label.grid(row=i, column=0, sticky="nw", padx=(0, 10), pady=1)
                            row_var = tk.StringVar(value="")
                            value_label = ttk.Label(
                                iono_box,
                                textvariable=row_var,
                                foreground=PALETTE["accent_dark"] if emphasis else PALETTE["text"],
                                font=(recap_font_family, 9, "bold") if emphasis else (recap_font_family, 9),
                                wraplength=620,
                                justify="left",
                            )
                            value_label.grid(row=i, column=1, sticky="w", pady=1)
                            widgets = [name_label, value_label]
                            if row_key == "guidance":
                                icon = _info_icon(iono_box, calculations.IONOSPHERIC_ZONES_SOURCE)
                                icon.grid(row=i, column=2, sticky="w", padx=(4, 0))
                                widgets.append(icon)
                            iono_rows[row_key] = (widgets, row_var)
                    elif kind == "eirp_mode":
                        mode_row = ttk.Frame(group_frame)
                        mode_row.grid(row=r, column=0, columnspan=4, sticky="w", pady=(0, 6))
                        ttk.Label(mode_row, text="EIRP source:", font=(recap_font_family, 9, "bold")).pack(
                            side=tk.LEFT, padx=(0, 8)
                        )
                        for value, text in (
                            ("eirp", f"EIRP value ({eirp_label.split(' [')[0]})"),
                            ("tx_chain", "Tx chain (power − losses + antenna gain)"),
                        ):
                            ttk.Radiobutton(
                                mode_row,
                                text=text,
                                value=value,
                                variable=eirp_mode_var,
                                command=lambda: _refresh_states(),
                            ).pack(side=tk.LEFT, padx=(0, 12))
                        _info_icon(
                            mode_row,
                            "Choose where the EIRP comes from: the EIRP value (shared with the main "
                            "window), or the Tx chain below (transmitter power − circuit loss − VSWR "
                            "loss + antenna gain). The fields of the other option are greyed out "
                            "and ignored; the EIRP actually used is shown under the Tx chain.",
                        ).pack(side=tk.LEFT)
                    elif kind == "gain_source":
                        gain_box = ttk.Frame(group_frame)
                        gain_box.grid(row=r, column=0, columnspan=4, sticky="w", pady=(2, 2))
                        gain_row = ttk.Frame(gain_box)
                        gain_row.pack(anchor="w")
                        ttk.Checkbutton(
                            gain_row,
                            text=(
                                "Antenna gain from the satellite antenna pattern at this elevation"
                                if rx_is_ground_station
                                else "G/T SAT at the antenna boresight: subtract the pattern roll-off "
                                "at this elevation"
                            ),
                            variable=gain_pattern_var,
                            command=lambda: _refresh_states(),
                        ).pack(side=tk.LEFT)
                        _info_icon(
                            gain_row,
                            (
                                "The satellite antenna pattern (gain versus angle from nadir, loaded "
                                "in the Orbit & Passes tab) gives the Tx antenna gain towards the "
                                "ground station: the angle from nadir follows from the elevation and "
                                "the satellite altitude. The pass analysis uses the same pattern at "
                                "every time step. The pattern already includes the off-pointing, so "
                                "the Tx Antenna Pointing Loss is not used."
                            )
                            if rx_is_ground_station
                            else (
                                "G/T SAT is usually specified at the antenna boresight (nadir), "
                                "where the satellite antenna gain is highest. Towards a station at "
                                "low elevation the antenna sees it far from nadir and its gain is "
                                "lower by the pattern roll-off (peak gain − gain at that angle), and "
                                "so is the G/T. With this option the G/T used is G/T SAT minus that "
                                "roll-off; the angle from nadir follows from the elevation and the "
                                "satellite altitude. The pass analysis does the same at every time "
                                "step. The pattern already includes the off-pointing, so the Rx "
                                "Antenna Depointing Angle is not used."
                            ),
                        ).pack(side=tk.LEFT, padx=(4, 8))
                        # Angle, gain and value used, on its own line under the option.
                        ttk.Label(
                            gain_box, textvariable=gain_pattern_text_var, foreground=PALETTE["accent_dark"]
                        ).pack(anchor="w", padx=(22, 0))
                    else:  # eirp_used
                        ttk.Label(
                            group_frame,
                            textvariable=eirp_used_var,
                            font=(recap_font_family, 9, "bold"),
                            foreground=PALETTE["accent_dark"],
                        ).grid(row=r, column=0, columnspan=4, sticky="w", pady=(4, 4))
                    r, c = r + 1, 0
                    continue
                padx = (0 if c == 0 else 18, 4)
                if kind == "scintillation":
                    scint_cell = ttk.Frame(group_frame)
                    scint_cell.grid(row=r, column=c * 2, columnspan=2, sticky="w", padx=padx, pady=2)
                    ttk.Checkbutton(
                        scint_cell,
                        text="Tropospheric scintillation",
                        variable=scint_var,
                    ).pack(side=tk.LEFT)
                    _info_icon(
                        scint_cell,
                        "Tropospheric scintillation fade (ITU-R P.618), added to the Rain Faded "
                        "column. The P.618 method is given for elevations of at least 5° and about "
                        "4–20 GHz: below that (e.g. L/S band at low elevation) it is an "
                        "extrapolation, and ionospheric scintillation dominates instead. Saved with "
                        "the inputs of this tab, independent of the switch of the pass analysis.",
                    ).pack(side=tk.LEFT, padx=(4, 0))
                    r, c = (r, 1) if c == 0 else (r + 1, 0)
                    continue
                if kind == "station":
                    ttk.Label(group_frame, text="Ground Station:").grid(
                        row=r, column=c * 2, sticky="w", padx=padx, pady=2
                    )
                    gs_cell = ttk.Frame(group_frame)
                    gs_cell.grid(row=r, column=c * 2 + 1, sticky="w", pady=2)
                    gs_box = ttk.Combobox(
                        gs_cell,
                        textvariable=station_var,
                        values=list(self.gs_menu["values"]) if self.gs_menu is not None else [],
                        state="readonly",
                        width=15,
                    )
                    gs_box.configure(
                        postcommand=lambda box=gs_box: box.configure(
                            values=list(self.gs_menu["values"]) if self.gs_menu is not None else []
                        )
                    )
                    gs_box.pack(side=tk.LEFT)
                    _info_icon(
                        gs_cell,
                        "Its location (latitude, longitude, altitude) sets the ITU-R atmospheric "
                        "attenuation; the geometry only uses elevation and altitude. Saved with "
                        "the inputs of this tab and independent of the station selected in the "
                        "Orbit & Passes tab (used as the default).",
                    ).pack(side=tk.LEFT, padx=(4, 0))
                    r, c = (r, 1) if c == 0 else (r + 1, 0)
                    continue
                if kind == "computed":
                    _kind, box_key, label_text = item
                    ttk.Label(group_frame, text=f"{label_text}:").grid(
                        row=r, column=c * 2, sticky="w", padx=padx, pady=2
                    )
                    cell = ttk.Frame(group_frame)
                    cell.grid(row=r, column=c * 2 + 1, sticky="w", pady=2)
                    box_var = tk.StringVar(value="")
                    ttk.Entry(cell, textvariable=box_var, width=15, state="readonly").pack(side=tk.LEFT)
                    box_icon = _info_icon(cell, "Computed from the coded bit rate and the overheads.")
                    box_icon.pack(side=tk.LEFT, padx=(4, 0))
                    computed_boxes[box_key] = (box_var, box_icon)
                    r, c = (r, 1) if c == 0 else (r + 1, 0)
                    continue
                if kind == "panel":
                    _kind, label_text, panel_entry, tooltip = item
                    label_widget = ttk.Label(group_frame, text=f"{label_text}:", foreground=PALETTE["accent"])
                    label_widget.grid(row=r, column=c * 2, sticky="w", padx=padx, pady=2)
                    cell = ttk.Frame(group_frame)
                    cell.grid(row=r, column=c * 2 + 1, sticky="w", pady=2)
                    linked = ttk.Entry(cell, textvariable=self._linked_var(panel_entry), width=15)
                    linked.pack(side=tk.LEFT)
                    linked.bind(
                        "<KeyRelease>", lambda event, edited=panel_entry: self.on_panel_entry_edited(edited)
                    )
                    panel_info = f"{tooltip} {shared_note}".strip()
                    panel_icon = _info_icon(cell, panel_info)
                    panel_icon.pack(side=tk.LEFT, padx=(4, 0))
                    info_icons[f"panel:{label_text}"] = (panel_icon, panel_info)
                    panel_specs.append((label_text, panel_entry))
                    field_widgets[f"panel:{label_text}"] = (label_widget, linked, PALETTE["accent"])
                    r, c = (r, 1) if c == 0 else (r + 1, 0)
                    continue
                _kind, key, label_text, default, placeholder, tooltip = item
                label_widget = ttk.Label(group_frame, text=f"{label_text}:")
                label_widget.grid(row=r, column=c * 2, sticky="w", padx=padx, pady=2)
                cell = ttk.Frame(group_frame)
                cell.grid(row=r, column=c * 2 + 1, sticky="w", pady=2)
                if placeholder:
                    entry = _PlaceholderEntry(cell, placeholder, width=15)
                else:
                    entry = ttk.Entry(cell, width=15)
                entry.pack(side=tk.LEFT)
                field_icon = _info_icon(cell, tooltip)
                field_icon.pack(side=tk.LEFT, padx=(4, 0))
                info_icons[key] = (field_icon, tooltip)
                entries[key] = entry
                entry_labels[key] = label_text
                entry_defaults[key] = default
                field_widgets[key] = (label_widget, entry, PALETTE["text"])
                r, c = (r, 1) if c == 0 else (r + 1, 0)

        tx_chain_keys = ("tx_power_w", "ant_circuit_loss_db", "vswr", "ant_gain_dbi")

        def _set_field_enabled(key: str, enabled: bool) -> None:
            label_widget, entry_widget, colour = field_widgets[key]
            entry_widget.state(["!disabled"] if enabled else ["disabled"])
            label_widget.configure(foreground=colour if enabled else PALETTE["text_muted"])
            # Placeholder entries set their own text colour: keep it in step.
            if isinstance(entry_widget, _PlaceholderEntry) and not entry_widget._showing_placeholder:
                entry_widget.configure(foreground=PALETTE["text"] if enabled else PALETTE["text_muted"])

        def _field_text(key: str) -> str:
            return field_widgets[key][1].get().strip()

        def _refresh_ionosphere() -> None:
            from_s4 = iono_mode_var.get() == "s4"
            _set_field_enabled("ionospheric_loss_db", not from_s4)
            _set_field_enabled("iono_s4_ref", from_s4)
            _set_field_enabled("iono_s4_ref_freq_mhz", from_s4)
            if not iono_rows:
                return
            values = {key: "" for key in iono_rows}
            station = station_var.get()
            if station in calculations.GROUND_STATIONS:
                lat, lon, _alt = calculations.GROUND_STATIONS[station]
                geomag = calculations.geomagnetic_latitude_deg(lat, lon)
                zone = calculations.ionospheric_zone(geomag)
                values["zone"] = f"{station} — geomagnetic latitude {geomag:.1f}° → {zone['zone']} (approx.)"
                values["guidance"] = zone["guidance"]
            if not from_s4:
                text = _field_text("ionospheric_loss_db") or "0"
                values["loss"] = f"{text} dB (manual value)"
            else:
                try:
                    s4_ref = float(_field_text("iono_s4_ref"))
                    ref_freq = float(_field_text("iono_s4_ref_freq_mhz"))
                    freq_mhz = float(_field_text(f"panel:{freq_label}")) * 1000.0
                    scint = calculations.ionospheric_scintillation_fade_db(s4_ref, freq_mhz, ref_freq)
                except ValueError:
                    values["loss"] = "enter S4, its reference frequency and the link frequency"
                else:
                    values["s4"] = (
                        f"{scint['s4']:.3f} at {freq_mhz:g} MHz "
                        f"(from {s4_ref:g} at {ref_freq:g} MHz, scaled as f^-1.5)"
                        + ("" if scint["valid"] else " — S4 ≥ 0.4: scaling not valid, indicative only")
                    )
                    values["fluctuation"] = f"{scint['peak_to_peak_db']:.2f} dB (27.5 · S4^1.26)"
                    values["loss"] = f"{scint['fade_db']:.2f} dB (half of the peak-to-peak fluctuation)"
            for key, (widgets, row_var) in iono_rows.items():
                row_var.set(values[key])
                # S4 rows only make sense when the loss is derived from S4.
                visible = from_s4 or key not in ("s4", "fluctuation")
                for widget in widgets:
                    if visible:
                        widget.grid()
                    else:
                        widget.grid_remove()

        def _ccsds_framing() -> tuple[float, dict]:
            """(coded bit rate [kbit/s], CCSDS overheads); raises ValueError."""

            try:
                coded_kbps = float(_field_text(f"panel:{bitrate_label}")) * 1000.0
                depth = float(_field_text("ccsds_interleaving"))
            except ValueError:
                raise ValueError(f"enter the {bitrate_label.split(' [')[0]} and the RS interleaving depth")
            if coded_kbps <= 0:
                raise ValueError("the coded bit rate must be greater than 0")
            packet_text = _field_text("ccsds_packet_length")
            try:
                packet_length = float(packet_text) if packet_text else None
                packet_sec_header = float(_field_text("ccsds_packet_sec_header") or 0)
            except ValueError:
                raise ValueError("the packet length and secondary header must be numbers")
            framing = calculations.ccsds_tm_overheads(
                depth,
                ccsds_ocf_var.get(),
                ccsds_fecf_var.get(),
                packet_length_bytes=packet_length,
                packet_secondary_header_bytes=packet_sec_header,
                packet_error_control_bytes=2 if ccsds_pec_var.get() else 0,
            )
            return coded_kbps, framing

        def _write_computed(key: str, text: str) -> None:
            """Fill a computed field, keeping it greyed out."""

            label_widget, entry_widget, _colour = field_widgets[key]
            if entry_widget.get().strip() == text:
                return
            entry_widget.state(["!disabled"])
            _set_entry_text(entry_widget, text)
            entry_widget.state(["disabled"])

        def _explain_framing(coded_kbps: float, framing: dict) -> None:
            """Put the breakdown of the CCSDS-derived values in their info icons."""

            depth = framing["interleaving_depth"]
            frame = framing["frame_bytes"]
            trailer = []
            trailer.append("OCF 4 B" if ccsds_ocf_var.get() else "no OCF")
            trailer.append("FECF 2 B" if ccsds_fecf_var.get() else "no FECF")
            packet_text = _field_text("ccsds_packet_length")
            if packet_text:
                packet_length = float(packet_text)
                cuc = float(_field_text("ccsds_packet_sec_header") or 0)
                pec = 2 if ccsds_pec_var.get() else 0
                per_packet = 6 + cuc + pec
                packet_line = (
                    f"  − packet overhead: header 6 + data field header (CUC) {cuc:g} + error "
                    f"control {pec} = {per_packet:g} B in every {packet_length:g}-B packet "
                    f"({framing['data_field_bytes'] - framing['user_data_bytes']:.1f} B per frame)\n"
                )
            else:
                packet_line = "  (packet overheads not deducted: no packet length given)\n"
            formatting_text = (
                f"Formatting factor = frame / user data = {frame} / {framing['user_data_bytes']:.1f} = "
                f"{framing['formatting']:.5f}\n"
                f"  Frame {frame} B (= 223 × I, I = {depth})\n"
                f"  − primary header 6 B (includes the First Header Pointer)\n"
                f"  − trailer: {', '.join(trailer)}\n"
                f"  = data field {framing['data_field_bytes']} B\n"
                f"{packet_line}"
                f"  = user data {framing['user_data_bytes']:.1f} B\n"
                "Computed from the CCSDS framing (Data framing above)."
            )
            coding_text = (
                f"Overhead = formatting × coding = {framing['formatting']:.5f} × {framing['coding']:.5f} "
                f"= {framing['total']:.5f}\n"
                f"  Coding = (RS codeblock {framing['codeblock_bytes']} B + ASM 4 B) / frame {frame} B × 2\n"
                f"  RS codeblock = frame {frame} B + {framing['codeblock_bytes'] - frame} B check symbols "
                f"(32 × I); CADU = {framing['codeblock_bytes'] + 4} B; × 2 = convolutional rate 1/2\n"
                f"Computed from the CCSDS framing; also used by the {passes_tab} tab."
            )
            bitrate_text = (
                f"Coded (channel) bit rate of the transmitter: {coded_kbps:.1f} kbit/s.\n"
                f"  Transfer frames = coded / coding = {coded_kbps:.1f} / {framing['coding']:.5f} = "
                f"{coded_kbps / framing['coding']:.1f} kbit/s (Eb/No reference: ASM, RS check "
                "symbols and convolutional coding removed)\n"
                f"  Information (user data) = coded / overhead = {coded_kbps:.1f} / {framing['total']:.5f} = "
                f"{coded_kbps / framing['total']:.1f} kbit/s\n"
                f"{shared_note}"
            )
            for key, text in (
                ("formatting_overhead", formatting_text),
                (f"panel:{overhead_label}", coding_text),
                (f"panel:{bitrate_label}", bitrate_text),
            ):
                info_icons[key][0].tooltip.text = text

        def _refresh_framing() -> None:
            if "ccsds_interleaving" not in field_widgets:
                return  # uplink: no CCSDS TM framing
            ccsds = framing_mode_var.get() == "ccsds"
            for key in ("ccsds_interleaving", "ccsds_packet_length", "ccsds_packet_sec_header"):
                _set_field_enabled(key, ccsds)
            for button in ccsds_option_buttons:
                button.state(["!disabled"] if ccsds else ["disabled"])
            for key in (f"panel:{overhead_label}", "formatting_overhead"):
                _set_field_enabled(key, not ccsds)
            explained = (f"panel:{bitrate_label}", f"panel:{overhead_label}", "formatting_overhead")

            def _restore_info() -> None:
                for key in explained:
                    icon, default_text = info_icons[key]
                    icon.tooltip.text = default_text

            if not ccsds:
                framing_used_var.set("")
                _restore_info()
                return
            try:
                coded_kbps, framing = _ccsds_framing()
            except ValueError as exc:
                framing_used_var.set(f"CCSDS framing: {exc}")
                _restore_info()
                return
            _explain_framing(coded_kbps, framing)
            _write_computed(f"panel:{overhead_label}", f"{framing['total']:.6g}")
            _write_computed("formatting_overhead", f"{framing['formatting']:.6g}")
            self.on_panel_entry_edited(overhead_entry)
            framing_used_var.set(
                f"Frame {framing['frame_bytes']} B (data field {framing['data_field_bytes']} B, user data "
                f"{framing['user_data_bytes']:.1f} B): "
                f"coding {framing['coding']:.5f} (ASM + RS(255,223) I={framing['interleaving_depth']} "
                f"+ conv 1/2), formatting {framing['formatting']:.5f}, total {framing['total']:.5f}"
            )

        def _refresh_rates() -> None:
            """Fill the computed rate boxes from the coded bit rate and overheads."""

            if not computed_boxes:
                return
            try:
                coded_kbps = float(_field_text(f"panel:{bitrate_label}")) * 1000.0
                overhead = float(_field_text(f"panel:{overhead_label}"))
                formatting = float(_field_text("formatting_overhead") or 1.0)
                if coded_kbps <= 0 or overhead <= 0 or formatting <= 0:
                    raise ValueError
            except ValueError:
                for box_var, box_icon in computed_boxes.values():
                    box_var.set("")
                    box_icon.tooltip.text = "Enter the coded bit rate, the overhead and the formatting factor."
                return
            coding = overhead / formatting
            frame_kbps = coded_kbps / coding
            info_kbps = coded_kbps / overhead
            values = {
                "coding_factor": (
                    f"{coding:.5f}",
                    f"Coding factor = overhead / formatting factor = {overhead:g} / {formatting:g} = "
                    f"{coding:.5f}: channel bits per transfer frame bit"
                    + (" (ASM + RS check symbols + convolutional 1/2)." if framing_mode_var.get() == "ccsds" else "."),
                ),
                "frame_rate": (
                    f"{frame_kbps:.1f}",
                    f"Transfer frame bit rate = coded / coding factor = {coded_kbps:.1f} / {coding:.5f} = "
                    f"{frame_kbps:.1f} kbit/s.\nThis is the bit rate Eb/No is referred to: the coded bit "
                    "rate with the channel coding (ASM, RS check symbols, convolutional code) removed.",
                ),
                "info_rate": (
                    f"{info_kbps:.1f}",
                    f"Information (user data) bit rate = coded / overhead = {coded_kbps:.1f} / {overhead:g} = "
                    f"{info_kbps:.1f} kbit/s\n= transfer frame bit rate / formatting factor: frame header and "
                    "trailer and packet overheads removed. It does not affect the margin.",
                ),
            }
            for key, (text, explanation) in values.items():
                box_var, box_icon = computed_boxes[key]
                box_var.set(text)
                box_icon.tooltip.text = explanation

        def _pattern_gain():
            """(angle from nadir, satellite antenna gain) at this view's geometry, or None."""

            try:
                elevation = float(_field_text("elevation_deg"))
                altitude = float(_field_text("sat_altitude_km"))
            except ValueError:
                return None
            station = calculations.GROUND_STATIONS.get(station_var.get())
            gs_alt_km = station[2] / 1000.0 if station else 0.0
            angle = calculations.off_nadir_angle_deg(elevation, altitude, gs_alt_km)
            return angle, calculations.satellite_antenna_gain_dbi(angle)

        def _pattern_rolloff():
            """(angle from nadir, pattern gain, peak gain − gain) at this view's geometry, or None."""

            angle_gain = _pattern_gain()
            if angle_gain is None:
                return None
            angle, gain = angle_gain
            return angle, gain, calculations.antenna_pattern(angle)

        def _refresh_states(_event=None) -> None:
            """Grey out the fields the calculation currently ignores."""

            _refresh_ionosphere()
            _refresh_framing()
            _refresh_rates()

            chain = eirp_mode_var.get() == "tx_chain"
            _set_field_enabled(f"panel:{eirp_label}", not chain)
            for key in tx_chain_keys:
                _set_field_enabled(key, chain)
            if rx_is_ground_station:
                # Tx antenna gain read from the satellite antenna pattern.
                use_pattern = bool(gain_pattern_var.get())
                gain_pattern_text_var.set("")
                if use_pattern:
                    _set_field_enabled("ant_gain_dbi", False)
                    angle_gain = _pattern_gain()
                    if angle_gain is None:
                        gain_pattern_text_var.set("enter Elevation Angle and Satellite Altitude")
                    else:
                        angle, gain = angle_gain
                        _write_computed("ant_gain_dbi", f"{gain:.2f}")
                        gain_pattern_text_var.set(f"{angle:.1f}° from nadir → {gain:.2f} dBi")
                _set_field_enabled("tx_pointing_loss_db", not use_pattern)
            else:
                # G/T SAT reduced by the satellite antenna pattern roll-off.
                gain_pattern_text_var.set("")
                use_pattern = bool(gain_pattern_var.get())
                if use_pattern:
                    rolloff = _pattern_rolloff()
                    if rolloff is None:
                        gain_pattern_text_var.set("enter Elevation Angle and Satellite Altitude")
                    else:
                        angle, gain, loss = rolloff
                        try:
                            gt_used = f"{float(_field_text(f'panel:{gt_label}')) - loss:.2f} dB/K"
                        except ValueError:
                            gt_used = "enter G/T SAT"
                        gain_pattern_text_var.set(
                            f"{angle:.1f}° from nadir → {gain:.2f} dBi, roll-off {loss:.2f} dB "
                            f"→ G/T used {gt_used}"
                        )
                _set_field_enabled("rx_depointing_deg", not use_pattern)
            # Axial ratios only matter when the polarisation loss is derived from them.
            derive_polarisation = not _field_text("polarisation_loss_db")
            _set_field_enabled("ant_axial_ratio_db", derive_polarisation)
            _set_field_enabled("rx_axial_ratio_db", derive_polarisation)
            _set_field_enabled("rx_beamwidth_deg", bool(_field_text("rx_depointing_deg")))

            if not chain:
                text = _field_text(f"panel:{eirp_label}")
                eirp_used_var.set(
                    f"EIRP used in the calculation: {text} dBW (EIRP value)"
                    if text
                    else f"EIRP used in the calculation: enter the {eirp_label.split(' [')[0]} value"
                )
                return
            try:
                power_w = float(_field_text("tx_power_w"))
                gain = float(_field_text("ant_gain_dbi"))
                circuit = float(_field_text("ant_circuit_loss_db") or 0.0)
                vswr_text = _field_text("vswr")
                vswr_loss = calculations.vswr_mismatch_loss_db(float(vswr_text)) if vswr_text else 0.0
                if power_w <= 0:
                    raise ValueError
            except ValueError:
                eirp_used_var.set(
                    "EIRP used in the calculation: enter Transmitter Power and Antenna Gain (Tx chain)"
                )
                return
            power_dbw = 10 * math.log10(power_w)
            eirp_used_var.set(
                f"EIRP used in the calculation: {power_dbw - circuit - vswr_loss + gain:.2f} dBW "
                f"(Tx chain: {power_dbw:.2f} dBW − {circuit:g} dB − {vswr_loss:.2f} dB VSWR "
                f"+ {gain:g} dBi)"
            )

        for key in tx_chain_keys + (
            "elevation_deg",
            "sat_altitude_km",
            "polarisation_loss_db",
            "rx_depointing_deg",
            f"panel:{eirp_label}",
            "ionospheric_loss_db",
            "iono_s4_ref",
            "iono_s4_ref_freq_mhz",
            f"panel:{freq_label}",
            f"panel:{bitrate_label}",
            f"panel:{overhead_label}",
            "formatting_overhead",
        ) + (
            ("ccsds_interleaving", "ccsds_packet_length", "ccsds_packet_sec_header")
            if rx_is_ground_station
            else ()
        ):
            field_widgets[key][1].bind("<KeyRelease>", _refresh_states, add="+")

        parameter_keys = {id(entry): key for key, entry in self._parameter_entries().items()}
        panel_keys = [
            (parameter_keys[id(panel_entry)], panel_entry)
            for _label, panel_entry in panel_specs
            if id(panel_entry) in parameter_keys
        ]

        tol_page = _add_page(
            "tolerances",
            "Tolerances",
            "Tolerances",
            "Optional favourable / adverse deviations from nominal, as positive dB values of effect on "
            "the margin. They are propagated through the link budget (Results › With Tolerances) and give "
            "the MEAN − 3σ and Worst Case RSS margins (requires a Required Eb/No).",
        )
        tol_frame = ttk.Frame(tol_page)
        tol_frame.pack(fill=tk.X)
        tolerance_help = {
            "Favourable [dB]": (
                "How much the margin can improve over nominal because of this contributor, as a "
                "positive dB value (e.g. +1.5 for an antenna gain of 0 ± 1.5 dBi)."
            ),
            "Adverse [dB]": (
                "How much the margin can get worse than nominal because of this contributor, as a "
                "positive dB value (e.g. 1.5 for an antenna gain of 0 ± 1.5 dBi)."
            ),
            "PDF": (
                "Probability distribution of the contributor between its adverse and favourable "
                "limits. UNI: any value equally likely (only the limits are known; most "
                "conservative). TRI: nominal most likely, limits rare (hardware with a measured "
                "nominal and specified limits). GAU: Gaussian with the limits at ±3σ (many small "
                "independent effects, e.g. propagation). It sets the mean and variance used for "
                "MEAN − 3σ; it does not affect Worst Case RSS."
            ),
        }
        for col, heading in enumerate(("Parameter", "Favourable [dB]", "Adverse [dB]", "PDF")):
            heading_cell = ttk.Frame(tol_frame)
            heading_cell.grid(row=1, column=col, sticky="w", padx=(0, 10), pady=1)
            ttk.Label(heading_cell, text=heading, font=(recap_font_family, 9, "bold")).pack(side=tk.LEFT)
            if heading in tolerance_help:
                _info_icon(heading_cell, tolerance_help[heading]).pack(side=tk.LEFT, padx=(4, 0))
        # (key, label, default distribution) -- defaults follow the reference
        # static link budget sheet (TRI for hardware, GAU for propagation).
        tolerance_specs = [
            ("tx_power", "Transmitter Power", "TRI"),
            ("tx_gain", "Tx Antenna Gain", "TRI"),
            ("propagation", "Atmospheric / Propagation Loss", "GAU"),
            ("rx_gt", "Rx G/T", "TRI"),
            ("rx_degradation", "Receiver Degradation", "TRI"),
        ]
        tolerance_widgets: dict[str, tuple[ttk.Entry, ttk.Entry, ttk.Combobox, str, str]] = {}
        for i, (key, label_text, default_pdf) in enumerate(tolerance_specs, start=2):
            ttk.Label(tol_frame, text=f"{label_text}:").grid(row=i, column=0, sticky="w", padx=(0, 10), pady=1)
            fav_entry = ttk.Entry(tol_frame, width=10)
            fav_entry.grid(row=i, column=1, sticky="w", padx=(0, 10), pady=1)
            adv_entry = ttk.Entry(tol_frame, width=10)
            adv_entry.grid(row=i, column=2, sticky="w", padx=(0, 10), pady=1)
            pdf_box = ttk.Combobox(
                tol_frame, values=list(calculations.TOLERANCE_DISTRIBUTIONS), width=6, state="readonly"
            )
            pdf_box.grid(row=i, column=3, sticky="w", pady=1)
            _add_tooltip(
                pdf_box,
                "TRI = triangular peaked at nominal, UNI = uniform between the limits, "
                "GAU = Gaussian with the limits taken as ±3σ.",
            )
            tolerance_widgets[key] = (fav_entry, adv_entry, pdf_box, label_text, default_pdf)

        def _format_value(value) -> str:
            if isinstance(value, bool) or value is None:
                return ""
            if isinstance(value, (int, float)):
                return format(float(value), ".10g")
            return str(value).strip()

        def _collect_values() -> dict:
            """Current inputs as a JSON-friendly dict (empty fields omitted)."""

            values: dict = {}
            for key, entry in entries.items():
                text = entry.get().strip()
                if not text:
                    continue
                try:
                    values[key] = float(text)
                except ValueError:
                    values[key] = text
            tolerances_out: dict = {}
            for key, (fav_entry, adv_entry, pdf_box, _label, default_pdf) in tolerance_widgets.items():
                item: dict = {}
                for name, entry in (("favourable_db", fav_entry), ("adverse_db", adv_entry)):
                    text = entry.get().strip()
                    if text:
                        try:
                            item[name] = float(text)
                        except ValueError:
                            item[name] = text
                if item or pdf_box.get() != default_pdf:
                    item["distribution"] = pdf_box.get()
                    tolerances_out[key] = item
            if tolerances_out:
                values["tolerances"] = tolerances_out
            values["eirp_mode"] = eirp_mode_var.get()
            values["ant_gain_from_pattern" if rx_is_ground_station else "gt_from_pattern"] = bool(
                gain_pattern_var.get()
            )
            values["station"] = station_var.get()
            values["ionospheric_mode"] = iono_mode_var.get()
            if rx_is_ground_station:
                values["framing_mode"] = framing_mode_var.get()
                values["ccsds_ocf"] = bool(ccsds_ocf_var.get())
                values["ccsds_fecf"] = bool(ccsds_fecf_var.get())
                values["ccsds_packet_error_control"] = bool(ccsds_pec_var.get())
            values["scintillation"] = bool(scint_var.get())
            return values

        def _apply_values(values: dict) -> None:
            """Fill the inputs from ``values``; anything missing goes back to its default."""

            for key in field_widgets:
                _set_field_enabled(key, True)  # disabled entries ignore edits
            for key, entry in entries.items():
                text = _format_value(values[key]) if key in values else entry_defaults[key]
                _set_entry_text(entry, text)
            mode = values.get("eirp_mode")
            if mode not in ("eirp", "tx_chain"):
                # Files saved before the explicit choice: the Tx chain was used
                # whenever a transmitter power was given.
                mode = "tx_chain" if _format_value(values.get("tx_power_w")) else "eirp"
            eirp_mode_var.set(mode)
            gain_pattern_var.set(
                bool(values.get("ant_gain_from_pattern" if rx_is_ground_station else "gt_from_pattern", False))
            )
            station = values.get("station")
            station_var.set(station if station in calculations.GROUND_STATIONS else self.gs_var.get())
            framing_mode = values.get("framing_mode")
            framing_mode_var.set(framing_mode if framing_mode in ("manual", "ccsds") else "manual")
            ccsds_ocf_var.set(bool(values.get("ccsds_ocf", True)))
            ccsds_fecf_var.set(bool(values.get("ccsds_fecf", False)))
            ccsds_pec_var.set(bool(values.get("ccsds_packet_error_control", True)))
            iono_mode = values.get("ionospheric_mode")
            iono_mode_var.set(iono_mode if iono_mode in ("manual", "s4") else "manual")
            scintillation = values.get("scintillation")
            scint_var.set(bool(scintillation) if isinstance(scintillation, bool) else self.scint_var.get())
            _refresh_states()
            _refresh_states()
            saved_tolerances = values.get("tolerances")
            if not isinstance(saved_tolerances, dict):
                saved_tolerances = {}
            for key, (fav_entry, adv_entry, pdf_box, _label, default_pdf) in tolerance_widgets.items():
                item = saved_tolerances.get(key)
                if not isinstance(item, dict):
                    item = {}
                _set_entry_text(fav_entry, _format_value(item.get("favourable_db")))
                _set_entry_text(adv_entry, _format_value(item.get("adverse_db")))
                distribution = str(item.get("distribution", default_pdf)).upper()
                pdf_box.set(distribution if distribution in calculations.TOLERANCE_DISTRIBUTIONS else default_pdf)
            _mark_stale()

        def _remember_values() -> None:
            self.static_lb_values[section] = _collect_values()

        # Results shown are "fresh" until any input they depend on changes.
        results_state = {"fresh": False}

        def _mark_stale(*_args) -> None:
            if not results_state["fresh"] or not win.winfo_exists():
                return
            results_state["fresh"] = False
            calc_button.configure(style="Tool.Accent.TButton")
            export_button.state(["disabled"])
            pdf_button.state(["disabled"])
            stale_label.pack(anchor="w", before=status_label)
            verdict_title.configure(text="⚠ OUTDATED", fg=PALETTE["warning"])
            for key in result_pages:
                if pages[key]["enabled"]:
                    _set_page_badge(key, "!", PALETTE["warning"])

        _apply_values(self.static_lb_values.get(section, {}))

        # Any edit of an input marks the results as out of date: the view's
        # own fields, the shared panel fields (also when edited in the main
        # window), the station, the EIRP source and the tolerances.
        for _label, entry_widget, _colour in field_widgets.values():
            entry_widget.bind("<KeyRelease>", _mark_stale, add="+")
        for fav_entry, adv_entry, pdf_box, _label, _default in tolerance_widgets.values():
            fav_entry.bind("<KeyRelease>", _mark_stale, add="+")
            adv_entry.bind("<KeyRelease>", _mark_stale, add="+")
            pdf_box.bind("<<ComboboxSelected>>", _mark_stale, add="+")
        station_var.trace_add("write", lambda *_args: _refresh_states())
        watched_vars = [
            station_var, eirp_mode_var, iono_mode_var, scint_var, framing_mode_var, ccsds_ocf_var, ccsds_fecf_var,
            ccsds_pec_var, gain_pattern_var,
        ] + [
            self._linked_var(panel_entry) for _label, panel_entry in panel_specs
        ]
        for var in watched_vars:
            var.trace_add("write", _mark_stale)

        # A parameters file loaded from the Orbit & Passes tab (or at
        # startup) refills this view.
        def _apply_loaded_values(values: dict) -> None:
            _apply_values(values)
            _mark_stale()

        self._static_lb_appliers[section] = _apply_loaded_values

        _nav_caption("RESULTS")
        result_pages = {
            "summary": "Summary",
            "link_budget": "Link Budget",
            "pfd": "PFD Check",
            "tolerance_results": "With Tolerances",
        }
        if not rx_is_ground_station:
            del result_pages["pfd"]  # the PFD limits apply to space-to-Earth links only
        for key, nav_text in result_pages.items():
            _add_page(key, nav_text)
            _set_page_enabled(key, False)
        _show_page("geometry")

        # Bottom of the sidebar: the key result of the last calculation, and
        # what the blue labels mean.
        sidebar_footer = tk.Frame(sidebar, bg=PALETTE["surface"])
        sidebar_footer.pack(side=tk.BOTTOM, fill=tk.X, padx=14, pady=14)
        verdict_card = tk.Frame(
            sidebar_footer, bg=PALETTE["surface_alt"], highlightthickness=1,
            highlightbackground=PALETTE["border"], padx=12, pady=10,
        )
        verdict_title = tk.Label(
            verdict_card, text="No results yet", bg=PALETTE["surface_alt"], fg=PALETTE["text_muted"],
            font=(recap_font_family, 9, "bold"), anchor="w", justify="left",
        )
        verdict_title.pack(fill=tk.X)
        verdict_value = tk.Label(
            verdict_card, text="Press Calculate", bg=PALETTE["surface_alt"], fg=PALETTE["text_muted"],
            font=(recap_font_family, 10), anchor="w", justify="left", wraplength=190,
        )
        verdict_value.pack(fill=tk.X, pady=(2, 0))
        verdict_card.pack(fill=tk.X)
        tk.Label(
            sidebar_footer,
            text=(
                "Blue labels: fields shared with the pass analysis tabs (editing them here "
                "changes them there too)."
            ),
            bg=PALETTE["surface"], fg=PALETTE["accent"], font=(recap_font_family, 8),
            anchor="w", justify="left", wraplength=200,
        ).pack(fill=tk.X, pady=(10, 0))
        license_text = self.license_status.summary()
        if license_text:
            tk.Label(
                sidebar_footer, text=license_text, bg=PALETTE["surface"],
                fg=PALETTE["warning"] if self.license_status.expiring_soon else PALETTE["text_muted"],
                font=(recap_font_family, 8), anchor="w", justify="left", wraplength=200,
            ).pack(fill=tk.X, pady=(8, 0))

        status_var = tk.StringVar(value="")
        stale_label = ttk.Label(
            status_box,
            text="⚠ Inputs changed: press Calculate to update the results",
            style="Toolbar.TLabel",
            foreground=PALETTE["danger"],
            font=(recap_font_family, 9, "bold"),
        )
        status_label = ttk.Label(
            status_box, textvariable=status_var, style="Toolbar.TLabel", foreground=PALETTE["text_muted"]
        )
        status_label.pack(anchor="w")
        status_box.bind(
            "<Configure>",
            lambda e: [label.configure(wraplength=max(e.width - 8, 120)) for label in (stale_label, status_label)],
        )

        def _set_verdict(checks) -> None:
            """Show the overall result in the sidebar card and the Summary entry."""

            if checks is None:
                verdict_title.configure(text="No results yet", fg=PALETTE["text_muted"])
                verdict_value.configure(text="Press Calculate", fg=PALETTE["text_muted"])
                return
            overall, text = static_report.summary_verdict(checks)
            margin = next((c for c in checks if c.label == "Data Recovery Margin [dB]"), None)
            colour = {"pass": PALETTE["success"], "fail": PALETTE["danger"]}.get(overall, PALETTE["text"])
            verdict_title.configure(
                text={"pass": "✓ PASS", "fail": "✗ FAIL"}.get(overall, "Not checked"), fg=colour
            )
            verdict_value.configure(
                text=(
                    f"Margin {margin.clear:.2f} dB clear / {margin.rain:.2f} dB rain"
                    if margin is not None
                    else text
                ),
                fg=PALETTE["text"],
            )
            if overall:
                _set_page_badge("summary", "✓" if overall == "pass" else "✗", colour)

        export_rows: list[tuple[str, str, str, str]] = []
        # Inputs and results of the last successful calculation, for the PDF report.
        last_report: dict = {}

        def _parse_required(entry: ttk.Entry, label: str) -> float:
            text = entry.get().strip()
            if not text:
                raise ValueError(f"'{label}' is empty. Please enter a value.")
            try:
                return float(text)
            except ValueError:
                raise ValueError(f"'{label}' must be a number (got '{text}').")

        def _parse_optional(entry: ttk.Entry, label: str, default: float = 0.0) -> float:
            text = entry.get().strip()
            if not text:
                return default
            try:
                return float(text)
            except ValueError:
                raise ValueError(f"'{label}' must be a number (got '{text}').")

        def _parse_optional_or_none(entry: ttk.Entry, label: str) -> float | None:
            text = entry.get().strip()
            if not text:
                return None
            try:
                return float(text)
            except ValueError:
                raise ValueError(f"'{label}' must be a number (got '{text}').")

        def _field(key: str, parser, **kwargs):
            return parser(entries[key], entry_labels[key], **kwargs)

        def _input_error(exc: Exception) -> None:
            """Show an input error on the page of the field it names, if any."""

            match = re.search(r"'([^']+)'", str(exc))
            widget = None
            if match:
                name = match.group(1)
                key = next((k for k, text in entry_labels.items() if text == name), f"panel:{name}")
                if key in field_widgets:
                    widget = field_widgets[key][1]
            page_key = _page_of(widget) if widget is not None else None
            if page_key:
                _show_page(page_key)
            messagebox.showerror("Input Error", str(exc), parent=win)
            if widget is not None and widget.winfo_exists():
                widget.focus_set()

        def _calculate():
            _remember_values()
            for key in result_pages:
                for widget in pages[key]["content"].winfo_children():
                    widget.destroy()
                _set_page_enabled(key, False)
                _set_page_badge(key)
            stale_label.pack_forget()
            _set_verdict(None)
            export_rows.clear()
            last_report.clear()
            export_button.state(["disabled"])
            pdf_button.state(["disabled"])
            try:
                elevation_deg = _field("elevation_deg", _parse_required)
                sat_altitude_km = _field("sat_altitude_km", _parse_required)
                pointing_loss_db = _field("tx_pointing_loss_db", _parse_optional)
                required_ebno = _field("required_ebno_db", _parse_optional_or_none)
                required_margin = _field("required_margin_db", _parse_optional_or_none)

                freq = _parse_required(freq_entry, freq_label) * u.GHz
                eirp = _parse_optional_or_none(eirp_entry, eirp_label)
                gt = _parse_required(gt_entry, gt_label)
                gt_pattern_rolloff_db = None
                if not rx_is_ground_station and gain_pattern_var.get():
                    rolloff = _pattern_rolloff()
                    if rolloff is None:
                        raise ValueError(
                            "'Elevation Angle [deg]' and 'Satellite Altitude [km]' are needed for "
                            "the G/T from the antenna pattern."
                        )
                    gt_pattern_rolloff_db = rolloff[2]
                    gt -= gt_pattern_rolloff_db
                demod_loss = _parse_required(demod_loss_entry, demod_loss_label)
                bitrate_mbps = _parse_required(bitrate_entry, bitrate_label)
                bitrate = bitrate_mbps * 1e6
                overhead = _parse_required(overhead_entry, overhead_label)
                other_att = _parse_optional(other_att_entry, other_att_label)
                link_availability = _parse_required(self.LA_entry, "Link Availability [%]")
                d_gs = _parse_required(self.d_gs_entry, "Antenna Diameter GS [m]")

                if eirp_mode_var.get() == "tx_chain":
                    tx_power_w = _field("tx_power_w", _parse_optional_or_none)
                    ant_circuit_loss_db = _field("ant_circuit_loss_db", _parse_optional)
                    vswr = _field("vswr", _parse_optional_or_none)
                    ant_gain_dbi = _field("ant_gain_dbi", _parse_optional_or_none)
                    if rx_is_ground_station and gain_pattern_var.get():
                        angle_gain = _pattern_gain()
                        if angle_gain is None:
                            raise ValueError(
                                "'Elevation Angle [deg]' and 'Satellite Altitude [km]' are needed for "
                                "the antenna gain from the pattern."
                            )
                        ant_gain_dbi = angle_gain[1]
                        pointing_loss_db = 0.0  # included in the pattern gain
                    if tx_power_w is None or ant_gain_dbi is None:
                        raise ValueError(
                            "The EIRP source is 'Tx chain': enter Transmitter Power [W] and "
                            "Antenna Gain [dBi], or choose 'EIRP value'."
                        )
                    if tx_power_w <= 0:
                        raise ValueError("'Transmitter Power [W]' must be greater than 0.")
                else:
                    # EIRP value: the Tx chain fields are ignored.
                    tx_power_w, ant_circuit_loss_db, vswr, ant_gain_dbi = None, 0.0, None, None
                    if eirp is None:
                        raise ValueError(
                            f"The EIRP source is 'EIRP value' but '{eirp_label}' is empty: enter it, "
                            "or choose 'Tx chain'."
                        )
                ant_axial_ratio_db = _field("ant_axial_ratio_db", _parse_optional_or_none)
                if iono_mode_var.get() == "s4":
                    ionospheric_loss_db = 0.0
                    iono_s4_ref = _field("iono_s4_ref", _parse_required)
                    iono_s4_ref_freq = _field("iono_s4_ref_freq_mhz", _parse_required)
                    if iono_s4_ref < 0 or iono_s4_ref_freq <= 0:
                        raise ValueError("S4 must be ≥ 0 and its reference frequency > 0.")
                else:
                    ionospheric_loss_db = _field("ionospheric_loss_db", _parse_optional)
                    iono_s4_ref, iono_s4_ref_freq = None, calculations.GPS_L1_MHZ
                polarisation_loss_db = _field("polarisation_loss_db", _parse_optional_or_none)
                rx_axial_ratio_db = _field("rx_axial_ratio_db", _parse_optional_or_none)
                multipath_loss_db = _field("multipath_loss_db", _parse_optional)
                modulation_degradation_db = _field("modulation_degradation_db", _parse_optional)
                pfd_limit = _field("pfd_limit", _parse_optional_or_none) if "pfd_limit" in entries else None
                ccsds_framing = None
                if rx_is_ground_station and framing_mode_var.get() == "ccsds":
                    # Validates the CCSDS inputs; the fields below were filled from them.
                    _coded_kbps, ccsds_framing = _ccsds_framing()
                formatting_overhead = _field("formatting_overhead", _parse_optional, default=1.0)
                if formatting_overhead <= 0:
                    raise ValueError("'Formatting Factor (frame / user data)' must be greater than 0.")
                rx_depointing_deg = _field("rx_depointing_deg", _parse_optional_or_none)
                if gt_pattern_rolloff_db is not None:
                    rx_depointing_deg = None  # included in the pattern roll-off
                rx_beamwidth_deg = _field("rx_beamwidth_deg", _parse_optional_or_none)
                if rx_beamwidth_deg is not None and rx_beamwidth_deg <= 0:
                    raise ValueError("'Rx Antenna 3dB Beamwidth [deg]' must be greater than 0.")
                if rx_depointing_deg is not None and rx_beamwidth_deg is None and not rx_is_ground_station:
                    raise ValueError(
                        "'Rx Antenna 3dB Beamwidth [deg]' is required when the Rx Antenna "
                        "Depointing Angle is set (it can only be derived from the GS antenna "
                        "diameter for the downlink)."
                    )

                tolerances = []
                for tol_key, (fav_entry, adv_entry, pdf_box, tol_label, _default_pdf) in tolerance_widgets.items():
                    fav = _parse_optional(fav_entry, f"{tol_label} favourable tolerance")
                    adv = _parse_optional(adv_entry, f"{tol_label} adverse tolerance")
                    if fav < 0 or adv < 0:
                        raise ValueError(
                            f"'{tol_label}' tolerances must be entered as positive dB values."
                        )
                    if fav or adv:
                        tolerances.append(
                            {
                                "name": tol_label,
                                "key": tol_key,
                                "favourable_db": fav,
                                "adverse_db": adv,
                                "distribution": pdf_box.get(),
                            }
                        )

                rolloff = _parse_optional_or_none(rolloff_entry, "Roll-off Factor")
                spectral_eff = _parse_optional_or_none(spectral_eff_entry, "Spectral Efficiency [bps/Hz]")
            except ValueError as exc:
                _input_error(exc)
                return

            gs_name = station_var.get()
            if not gs_name or gs_name not in calculations.GROUND_STATIONS:
                messagebox.showerror("Error", "Please select a valid Ground Station first.", parent=win)
                return
            lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]

            occupied_bandwidth_hz = None
            if rolloff is not None and spectral_eff:
                occupied_bandwidth_hz = bitrate_mbps * (1 + rolloff) / spectral_eff * 1e6

            info_bitrate_kbps = bitrate_mbps * 1000.0 / overhead if overhead else None

            try:
                budget = calculate_fixed_elevation_link_budget(
                    freq=freq,
                    elevation_deg=elevation_deg,
                    sat_altitude_km=sat_altitude_km,
                    lat_gs=lat_gs,
                    lon_gs=lon_gs,
                    alt_gs_km=alt_gs_m / 1000.0,
                    d_gs=d_gs,
                    eirp=eirp,
                    gt=gt,
                    demod_loss=demod_loss,
                    bitrate=bitrate,
                    overhead=overhead,
                    other_att=other_att,
                    pointing_loss_db=pointing_loss_db,
                    link_availability_pct=link_availability,
                    include_scintillation=scint_var.get(),
                    required_ebno=required_ebno,
                    ionospheric_s4_ref=iono_s4_ref,
                    ionospheric_s4_ref_freq_mhz=iono_s4_ref_freq,
                    tx_power_w=tx_power_w,
                    antenna_circuit_loss_db=ant_circuit_loss_db,
                    vswr=vswr,
                    antenna_gain_dbi=ant_gain_dbi,
                    ionospheric_loss_db=ionospheric_loss_db,
                    polarisation_loss_db=polarisation_loss_db,
                    multipath_loss_db=multipath_loss_db,
                    modulation_degradation_db=modulation_degradation_db,
                    occupied_bandwidth_hz=occupied_bandwidth_hz,
                    pfd_limit_dbw_m2_4khz=pfd_limit,
                    tx_axial_ratio_db=ant_axial_ratio_db,
                    rx_axial_ratio_db=rx_axial_ratio_db,
                    rx_depointing_deg=rx_depointing_deg,
                    rx_beamwidth_3db_deg=rx_beamwidth_deg,
                    formatting_overhead=formatting_overhead,
                    tolerances=tolerances or None,
                )
            except ValueError as exc:
                _input_error(exc)
                return
            extra_info = {
                "station": f"{gs_name} ({lat_gs:.4f}°, {lon_gs:.4f}°, {alt_gs_m:g} m)",
                "link_availability_pct": link_availability,
                "formatting_overhead": formatting_overhead,
                "info_bitrate_kbps": info_bitrate_kbps,
                "bitrate_kbps": bitrate_mbps * 1000.0,
                "code_rate": (1.0 / overhead) if overhead else None,
                "occupied_bandwidth_khz": occupied_bandwidth_hz / 1000.0 if occupied_bandwidth_hz else None,
                "ant_axial_ratio_db": ant_axial_ratio_db,
                "rx_axial_ratio_db": rx_axial_ratio_db,
                "pointing_loss_db": pointing_loss_db,
                "required_margin_db": required_margin,
                "bandwidth_missing": occupied_bandwidth_hz is None,
                "pfd_applicable": rx_is_ground_station,
                "elevation_deg": elevation_deg,
                "tolerances_ignored": bool(tolerances) and required_ebno is None,
                "ccsds_framing": ccsds_framing,
            }
            # Inputs actually used, in the same sections as this view (CSV and PDF).
            input_groups = []
            for group_title, items in field_groups:
                group_items = []
                for item in items:
                    if item[0] == "station":
                        group_items.append(("Ground Station", gs_name))
                    elif item[0] == "scintillation":
                        group_items.append(
                            ("Tropospheric Scintillation", "Included" if scint_var.get() else "Excluded")
                        )
                    elif item[0] == "iono_mode":
                        group_items.append(
                            (
                                "Ionospheric Loss",
                                "From S4 (ITU-R P.531)" if iono_mode_var.get() == "s4" else "Manual",
                            )
                        )
                    elif item[0] == "framing_mode":
                        group_items.append(
                            (
                                "Data Framing",
                                "CCSDS TM (RS(255,223) + conv 1/2)" if framing_mode_var.get() == "ccsds" else "Manual",
                            )
                        )
                    elif item[0] == "ccsds_options":
                        if framing_mode_var.get() == "ccsds":
                            group_items.append(
                                (
                                    "Frame OCF / FECF, Packet Error Control",
                                    f"{'yes' if ccsds_ocf_var.get() else 'no'} / "
                                    f"{'yes' if ccsds_fecf_var.get() else 'no'}, "
                                    f"{'yes' if ccsds_pec_var.get() else 'no'}",
                                )
                            )
                    elif item[0] == "gain_source" and not rx_is_ground_station:
                        if gt_pattern_rolloff_db is not None:
                            angle_gain = _pattern_gain()
                            group_items.append(
                                (
                                    "G/T SAT Pattern Roll-off [dB]",
                                    f"{gt_pattern_rolloff_db:.2f} ({angle_gain[0]:.1f}° from nadir; "
                                    f"G/T used {gt:.2f} dB/K)",
                                )
                            )
                    elif item[0] == "gain_source":
                        if eirp_mode_var.get() == "tx_chain" and gain_pattern_var.get():
                            angle_gain = _pattern_gain()
                            group_items.append(
                                (
                                    "Antenna Gain Source",
                                    f"Satellite antenna pattern ({angle_gain[0]:.1f}° from nadir)"
                                    if angle_gain
                                    else "Satellite antenna pattern",
                                )
                            )
                    elif item[0] in ("iono_used", "eirp_used", "framing_used", "computed"):
                        continue
                    elif item[0] == "eirp_mode":
                        group_items.append(
                            (
                                "EIRP Source",
                                "Tx chain" if eirp_mode_var.get() == "tx_chain" else "EIRP value",
                            )
                        )
                    elif item[0] in ("panel", "field"):
                        key = f"panel:{item[1]}" if item[0] == "panel" else item[1]
                        label_widget, entry_widget, _colour = field_widgets[key]
                        text = entry_widget.get().strip()
                        # Fields greyed out (ignored by the calculation) are left out.
                        if text and not entry_widget.instate(["disabled"]):
                            label_text = item[1] if item[0] == "panel" else entry_labels[key]
                            group_items.append((label_text, text))
                input_groups.append((group_title, group_items))
            project = self.project_var.get().strip()
            if project:
                export_rows.append(("Project", project, "", ""))
            for group_title, group_items in input_groups:
                for label_text, text in group_items:
                    export_rows.append((f"Input: {group_title}", label_text, text, ""))
            export_rows.extend(
                self._render_fixed_link_budget_results(
                    {key: pages[key]["content"] for key in result_pages},
                    budget,
                    link_availability,
                    extra_info,
                )
            )
            for key in result_pages:
                _set_page_enabled(key, bool(pages[key]["content"].winfo_children()))
            _set_verdict(static_report.summary_checks(budget, extra_info))
            export_button.state(["!disabled"])

            last_report.update(
                budget=budget,
                extra_info=extra_info,
                context=static_report.ReportContext(
                    title=title,
                    direction=direction_text,
                    station=gs_name,
                    elevation_deg=elevation_deg,
                    link_availability_pct=link_availability,
                    input_groups=input_groups,
                    logo_path=_asset_path(LOGO_FILE),
                    project=self.project_var.get().strip(),
                    app_version=APP_VERSION,
                ),
            )
            pdf_button.state(["!disabled"])
            results_state["fresh"] = True
            calc_button.configure(style="Tool.TButton")
            if current_page["key"] not in result_pages:
                _show_page("summary")

        def _default_parameters_path() -> str:
            return (
                self.current_param_file
                or calculations.resolve_optional_data_file("parameters.json", env_var="PARAMETERS_FILE")
                or os.path.join(os.getcwd(), "parameters.json")
            )

        def _save():
            default_path = _default_parameters_path()
            file_path = filedialog.asksaveasfilename(
                parent=win,
                title="Save Fixed-Elevation Parameters",
                initialdir=os.path.dirname(os.path.abspath(default_path)),
                initialfile=os.path.basename(default_path),
                defaultextension=".json",
                filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
            )
            if not file_path:
                return
            _remember_values()
            panel_values = {}
            for key, panel_entry in panel_keys:
                text = panel_entry.get().strip()
                if text:
                    try:
                        panel_values[key] = float(text)
                    except ValueError:
                        panel_values[key] = text
            project = self.project_var.get().strip()
            if project:
                panel_values["project_name"] = project
            norad = self.norad_var.get().strip()
            if norad.isdigit():
                panel_values["norad_id"] = int(norad)
            orbit = self._custom_orbit_for_file()
            if orbit is not None:
                panel_values["orbit"] = orbit
            try:
                calculations.save_json_section(
                    file_path, section, self.static_lb_values[section], top_level=panel_values
                )
            except Exception as exc:
                messagebox.showerror("Save Error", f"Unable to save the parameters: {exc}", parent=win)
                return
            status_var.set(
                f"Saved to {os.path.abspath(file_path)} (section '{section}', the shared values and the "
                "project name)"
            )

        def _load():
            file_path = filedialog.askopenfilename(
                parent=win,
                title="Load Fixed-Elevation Parameters",
                initialdir=os.path.dirname(os.path.abspath(_default_parameters_path())),
                filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
            )
            if not file_path:
                return
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    payload = json.load(f)
                values = calculations.extract_json_section(payload, section)
            except Exception as exc:
                messagebox.showerror("Load Error", f"Unable to load the parameters: {exc}", parent=win)
                return
            for key, panel_entry in panel_keys:
                if key in payload:
                    panel_entry.delete(0, tk.END)
                    panel_entry.insert(0, _format_value(payload[key]))
                    self.on_panel_entry_edited(panel_entry)
            _apply_values(values)
            _remember_values()
            status_var.set(f"Loaded from {os.path.abspath(file_path)} (section '{section}')")

        def _reset():
            _apply_values({})
            _remember_values()
            status_var.set("Inputs reset to their defaults.")

        def _export_csv():
            if not export_rows:
                return
            file_path = filedialog.asksaveasfilename(
                parent=win,
                title="Export Fixed-Elevation Link Budget",
                defaultextension=".csv",
                filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
            )
            if not file_path:
                return
            try:
                with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Section", "Parameter", "Clear Sky / Value", "Rain Faded"])
                    writer.writerows(export_rows)
            except Exception as exc:
                messagebox.showerror("Export Error", f"Failed to export CSV: {exc}", parent=win)
                return
            status_var.set(f"Results exported to {os.path.abspath(file_path)}")

        def _export_pdf():
            if not last_report:
                return
            direction_slug = "downlink" if rx_is_ground_station else "uplink"
            file_path = filedialog.asksaveasfilename(
                parent=win,
                title="Export Fixed-Elevation Link Budget Report",
                initialfile=(
                    _project_slug(self.project_var.get())
                    + f"static_link_budget_{direction_slug}_{datetime.now():%Y%m%d}.pdf"
                ),
                defaultextension=".pdf",
                filetypes=[("PDF Files", "*.pdf"), ("All Files", "*.*")],
            )
            if not file_path:
                return
            try:
                static_report.build_static_report_pdf(
                    file_path, last_report["budget"], last_report["extra_info"], last_report["context"]
                )
            except ImportError:
                messagebox.showerror(
                    "Export Error",
                    "The PDF report needs the 'reportlab' package (pip install reportlab).",
                    parent=win,
                )
                return
            except Exception as exc:
                messagebox.showerror("Export Error", f"Failed to create the PDF report: {exc}", parent=win)
                return
            status_var.set(f"PDF report saved to {os.path.abspath(file_path)}")

        # Toolstrip: groups of large buttons, each with its caption below, as
        # in the MATLAB toolstrip.
        def _tool_group(caption: str, first: bool = False) -> ttk.Frame:
            if not first:
                ttk.Separator(button_row, orient="vertical").pack(side=tk.LEFT, fill=tk.Y, padx=10)
            group = ttk.Frame(button_row, style="Toolbar.TFrame")
            group.pack(side=tk.LEFT, fill=tk.Y)
            buttons = ttk.Frame(group, style="Toolbar.TFrame")
            buttons.pack(side=tk.TOP)
            ttk.Label(
                group, text=caption, style="ToolCaption.TLabel"
            ).pack(side=tk.BOTTOM, pady=(4, 0))
            return buttons

        run_group = _tool_group("RUN", first=True)
        file_group = _tool_group("FILE")
        export_group = _tool_group("EXPORT")

        calc_button = ttk.Button(run_group, text="Calculate", command=_calculate, style="Tool.Accent.TButton")
        calc_button.pack(side=tk.LEFT)
        _add_tooltip(calc_button, "Blue when there are no results yet or the inputs changed since the last calculation.")
        load_button = ttk.Button(file_group, text="Load...", command=_load, style="Tool.TButton")
        load_button.pack(side=tk.LEFT, padx=(0, 6))
        _add_tooltip(
            load_button,
            f"Load this tab's inputs from the '{section}' section of a parameters JSON file, "
            "and the shared values (blue labels) when the file has them.",
        )
        save_button = ttk.Button(file_group, text="Save...", command=_save, style="Tool.TButton")
        save_button.pack(side=tk.LEFT, padx=(0, 6))
        _add_tooltip(
            save_button,
            f"Save this tab's inputs as the '{section}' section of a JSON file, together "
            "with the shared values (blue labels); the other keys already in "
            "the file are kept. Saving into parameters.json makes them load automatically at "
            "startup.",
        )
        ttk.Button(file_group, text="Reset", command=_reset, style="Tool.TButton").pack(side=tk.LEFT)
        export_button = ttk.Button(export_group, text="Export CSV...", command=_export_csv, style="Tool.TButton")
        export_button.pack(side=tk.LEFT, padx=(0, 6))
        export_button.state(["disabled"])
        pdf_button = ttk.Button(export_group, text="Export PDF...", command=_export_pdf, style="Tool.TButton")
        pdf_button.pack(side=tk.LEFT)
        pdf_button.state(["disabled"])
        _add_tooltip(
            pdf_button,
            "Create a PDF report (A4) of the last calculation: summary with pass/fail checks, "
            "inputs, link budget table, PFD check, tolerances and method.",
        )

        def _pass_link_losses() -> dict:
            """This tab's elevation-independent losses, for the pass analysis.

            Raises ValueError (naming the field) when an input is invalid.
            """

            def opt(key, parser, **kwargs):
                # Fields that only one direction has count as empty in the other.
                return _field(key, parser, **kwargs) if key in entries else kwargs.get("default")

            freq_ghz = _parse_required(freq_entry, freq_label)
            d_gs = _parse_required(self.d_gs_entry, "Antenna Diameter GS [m]")
            if iono_mode_var.get() == "s4":
                iono_kwargs = dict(
                    ionospheric_s4_ref=_field("iono_s4_ref", _parse_required),
                    ionospheric_s4_ref_freq_mhz=_field("iono_s4_ref_freq_mhz", _parse_required),
                )
            else:
                iono_kwargs = dict(ionospheric_loss_db=opt("ionospheric_loss_db", _parse_optional, default=0.0))
            losses = calculations.static_link_losses(
                freq_ghz,
                d_gs,
                polarisation_loss_db=opt("polarisation_loss_db", _parse_optional_or_none),
                tx_axial_ratio_db=opt("ant_axial_ratio_db", _parse_optional_or_none),
                rx_axial_ratio_db=opt("rx_axial_ratio_db", _parse_optional_or_none),
                multipath_loss_db=opt("multipath_loss_db", _parse_optional, default=0.0),
                rx_depointing_deg=opt("rx_depointing_deg", _parse_optional_or_none),
                rx_beamwidth_3db_deg=opt("rx_beamwidth_deg", _parse_optional_or_none),
                **iono_kwargs,
            )
            formatting = opt("formatting_overhead", _parse_optional, default=1.0)
            if formatting <= 0:
                raise ValueError("'Formatting Factor (frame / user data)' must be greater than 0.")
            tx_chain_dbw = None
            if eirp_mode_var.get() == "tx_chain":
                power_w = _field("tx_power_w", _parse_required)
                if power_w <= 0:
                    raise ValueError("'Transmitter Power [W]' must be greater than 0.")
                vswr = _field("vswr", _parse_optional_or_none)
                tx_chain_dbw = (
                    10 * math.log10(power_w)
                    - _field("ant_circuit_loss_db", _parse_optional)
                    - (calculations.vswr_mismatch_loss_db(vswr) if vswr else 0.0)
                )
            return {
                "tx_chain_dbw": tx_chain_dbw,
                "ionospheric_loss_db": losses["ionospheric_loss_db"],
                "polarisation_loss_db": losses["polarisation_loss_db"],
                "multipath_loss_db": losses["multipath_loss_db"],
                "rx_pointing_loss_db": losses["rx_pointing_loss_db"],
                "tx_pointing_loss_db": opt("tx_pointing_loss_db", _parse_optional, default=0.0),
                "modulation_degradation_db": opt("modulation_degradation_db", _parse_optional, default=0.0),
                "formatting_overhead": formatting,
            }

        self._static_lb_loss_providers[section] = _pass_link_losses
        self._static_lb_refreshers.append(_refresh_states)


    def _render_fixed_link_budget_results(self, pages, budget, link_availability_pct, extra_info=None):
        """Render the fixed-elevation results and return them as CSV export rows.

        ``pages`` maps the result pages of the static view ("summary",
        "link_budget", optionally "pfd", "tolerance_results") to the frames
        to fill: the summary (overall verdict, key results checked against
        their requirements, notes and the static-budget info), the
        two-condition (clear/rain) link table, the power flux density check
        with an explanation of how it is computed, and the link budget with
        tolerances. Pages with nothing to show are left empty. The rows come
        from :mod:`static_report`, shared with the PDF report. Every table
        row has a hover tooltip with its definition, and margins are marked
        pass/fail against their requirements. Returns
        ``[(section, parameter, clear_or_value, rain), ...]``.
        """

        extra_info = extra_info or {}
        export_rows: list[tuple[str, str, str, str]] = []
        rain_p = static_report.rain_percentage(link_availability_pct)
        font_family = _pick_font(["Segoe UI", "Helvetica Neue", "Helvetica", "Arial"])

        def page_heading(parent, heading: str, description: str = "") -> None:
            ttk.Label(parent, text=heading, style="PageTitle.TLabel").pack(anchor="w")
            if description:
                ttk.Label(
                    parent, text=description, style="PageSub.TLabel", wraplength=760, justify="left"
                ).pack(anchor="w", pady=(2, 0))
            ttk.Frame(parent, height=12).pack(fill=tk.X)

        # ---- Summary ----------------------------------------------------
        summary = pages["summary"]
        page_heading(
            summary,
            "Summary",
            f"Slant range {budget['slant_range_km']:.1f} km  ·  free space loss {budget['path_loss_db']:.2f} dB  ·  "
            f"rain faded at p = {rain_p:g}% of the time",
        )
        export_rows.append(("Geometry", "Slant Range [km]", f"{budget['slant_range_km']:.2f}", ""))

        checks = static_report.summary_checks(budget, extra_info)
        overall, verdict_text = static_report.summary_verdict(checks)
        banner_bg, banner_fg = {
            "pass": (PALETTE["success_light"], PALETTE["success"]),
            "fail": (PALETTE["danger_light"], PALETTE["danger"]),
        }.get(overall, (PALETTE["surface_alt"], PALETTE["text_muted"]))
        banner = tk.Frame(summary, bg=banner_bg, highlightthickness=1, highlightbackground=banner_fg)
        banner.pack(fill=tk.X)
        tk.Label(
            banner,
            text={"pass": "✓  PASS", "fail": "✗  FAIL"}.get(overall, "NOT CHECKED"),
            bg=banner_bg, fg=banner_fg, font=(font_family, 14, "bold"), padx=16, pady=10,
        ).pack(side=tk.LEFT)
        tk.Label(
            banner, text=verdict_text, bg=banner_bg, fg=PALETTE["text"], font=(font_family, 10),
            wraplength=560, justify="left", anchor="w",
        ).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 16))
        export_rows.append(("Summary", verdict_text, "", ""))

        if checks:
            tiles = ttk.Frame(summary)
            tiles.pack(fill=tk.X, pady=(14, 0))
            per_row = 3
            for col in range(per_row):
                tiles.grid_columnconfigure(col, weight=1, uniform="tile")
            for i, check in enumerate(checks):
                r, c = divmod(i, per_row)
                tile = tk.Frame(
                    tiles, bg=PALETTE["surface"], highlightthickness=1, highlightbackground=PALETTE["border"],
                    padx=14, pady=10,
                )
                tile.grid(row=r, column=c, sticky="nsew", padx=(0 if c == 0 else 10, 0), pady=(0, 10))
                status_colour = {"pass": PALETTE["success"], "fail": PALETTE["danger"]}.get(
                    check.status, PALETTE["text_muted"]
                )
                top = tk.Frame(tile, bg=PALETTE["surface"])
                top.pack(fill=tk.X)
                tk.Label(
                    top, text={"pass": "PASS", "fail": "FAIL"}.get(check.status, "–"), bg=status_colour
                    if check.status else PALETTE["surface"], fg="white" if check.status else PALETTE["text_muted"],
                    font=(font_family, 8, "bold"), padx=6,
                ).pack(side=tk.RIGHT, anchor="n")
                tk.Label(
                    top, text=check.label.replace("MARGIN – ", "Margin – "), bg=PALETTE["surface"],
                    fg=PALETTE["text_muted"], font=(font_family, 9, "bold"), anchor="w", justify="left",
                    wraplength=170,
                ).pack(side=tk.LEFT)
                same = check.clear == check.rain
                value_text = f"{check.clear:.2f}" if same else f"{check.clear:.2f}  /  {check.rain:.2f}"
                tk.Label(
                    tile, text=value_text, bg=PALETTE["surface"], fg=status_colour if check.status else PALETTE["text"],
                    font=(font_family, 17, "bold"), anchor="w",
                ).pack(fill=tk.X, pady=(4, 0))
                tk.Label(
                    tile,
                    text=("" if same else f"clear sky / rain faded  ·  ") + f"required {check.requirement}",
                    bg=PALETTE["surface"], fg=PALETTE["text_muted"], font=(font_family, 8), anchor="w",
                ).pack(fill=tk.X)

        def section_heading(parent, text: str) -> None:
            ttk.Label(parent, text=text, style="Section.TLabel").pack(anchor="w", pady=(14, 6))

        notes = static_report.analysis_notes(extra_info, budget)
        if notes:
            section_heading(summary, "Notes")
            notes_card = tk.Frame(
                summary, bg=PALETTE["danger_light"], highlightthickness=1,
                highlightbackground=PALETTE["danger_border"], padx=14, pady=8,
            )
            notes_card.pack(fill=tk.X)
            for note in notes:
                tk.Label(
                    notes_card, text=f"•  {note}", bg=PALETTE["danger_light"], fg=PALETTE["danger_dark"],
                    font=(font_family, 9), wraplength=740, justify="left", anchor="w",
                ).pack(fill=tk.X, pady=2)
                export_rows.append(("Note", note, "", ""))

        info_rows = static_report.info_rows(budget, extra_info)
        if info_rows:
            section_heading(summary, "Static Link Budget Info")
            info_table = ttk.Treeview(
                summary, columns=("param", "value"), show="headings", height=len(info_rows), selectmode="none"
            )
            info_table.heading("param", text="Parameter")
            info_table.heading("value", text="Value")
            info_table.column("param", anchor="w", width=380)
            info_table.column("value", anchor="w", width=420)
            _configure_result_tags(info_table)
            for i, (label_text, value_text) in enumerate(info_rows):
                info_table.insert(
                    "", "end", values=(label_text, value_text), tags=("evenrow" if i % 2 == 0 else "oddrow",)
                )
                export_rows.append(("Info", label_text, value_text, ""))
            info_table.pack(fill=tk.X)

        # ---- Link budget ------------------------------------------------
        parent = pages["link_budget"]
        page_heading(
            parent,
            "Link Budget",
            "Clear sky and rain faded conditions. Hover over a row for its definition. Margins are "
            "marked ✓ / ✗ against their requirement (nominal margin ≥ Required Nominal Margin, "
            "MEAN − 3σ and Worst Case RSS ≥ 0 dB).",
        )

        required_margin = extra_info.get("required_margin_db")
        rows = static_report.link_budget_rows(budget, required_margin)

        table_frame = ttk.Frame(parent)
        table_frame.pack(fill=tk.BOTH, expand=True)
        groups = list(dict.fromkeys(row.group for row in rows))
        table = ttk.Treeview(
            table_frame, columns=("param", "clear", "rain"), show="headings", height=len(rows) + len(groups)
        )
        table.heading("param", text="Parameter")
        table.heading("clear", text="Clear Sky")
        table.heading("rain", text=f"Rain Faded (p={rain_p:g}%)")
        table.column("param", anchor="w", width=360)
        table.column("clear", anchor="center", width=150)
        table.column("rain", anchor="center", width=190)
        _configure_result_tags(table)

        row_help: dict[str, str] = {}
        current_group = None
        i = 0
        for row in rows:
            if row.group != current_group:
                current_group = row.group
                table.insert("", "end", values=(row.group.upper(), "", ""), tags=("group",))
                i = 0
            clear_status = row.status(row.clear)
            rain_status = row.status(row.rain)
            clear_text = _format_checked(row.clear, clear_status)
            rain_text = _format_checked(row.rain, rain_status)
            tags = ["evenrow" if i % 2 == 0 else "oddrow"]
            i += 1
            if row.emphasis:
                tags.append("emphasis")
            if "fail" in (clear_status, rain_status):
                tags.append("fail")
            elif "pass" in (clear_status, rain_status):
                tags.append("pass")
            iid = table.insert("", "end", values=(row.label, clear_text, rain_text), tags=tuple(tags))
            help_text = _FIXED_LB_ROW_HELP.get(row.label)
            if row.label == "Data Recovery Margin [dB]" and required_margin is not None:
                help_text = f"{help_text} Requirement: ≥ {required_margin:g} dB."
            if help_text:
                row_help[iid] = help_text
            export_rows.append(("Link Budget", row.label, clear_text, rain_text))
        _TreeviewRowTooltip(table, row_help)
        table.pack(fill=tk.BOTH, expand=True)

        if "pfd" in pages and static_report.pfd_applicable(extra_info):
            page_heading(pages["pfd"], "Power Flux Density", static_report.PFD_EXPLANATION)
            self._render_pfd_section(pages["pfd"], budget, extra_info, export_rows)

        # Tolerances propagated through the link budget (worst case: rain faded).
        tol_rows = static_report.tolerance_rows(budget, "rain_faded")
        if tol_rows and "tolerance_results" in pages:
            tol_frame = pages["tolerance_results"]
            page_heading(
                tol_frame,
                "Link Budget with Tolerances",
                f"Rain faded condition (p = {rain_p:g}%), the worst case. "
                + static_report.TOLERANCE_TABLE_EXPLANATION,
            )
            columns = ("param", "nominal", "fav", "adv", "mean", "var", "pdf")
            tol_table = ttk.Treeview(tol_frame, columns=columns, show="headings", height=len(tol_rows))
            for col, heading, width, anchor in (
                ("param", "Parameter", 270, "w"),
                ("nominal", "Nominal", 90, "center"),
                ("fav", "Favourable", 105, "center"),
                ("adv", "Adverse", 95, "center"),
                ("mean", "Mean", 85, "center"),
                ("var", "Variance [dB²]", 125, "center"),
                ("pdf", "PDF", 90, "center"),
            ):
                tol_table.heading(col, text=heading)
                tol_table.column(col, anchor=anchor, width=width)
            _configure_result_tags(tol_table)
            for i, row in enumerate(tol_rows):
                tags = ("evenrow" if i % 2 == 0 else "oddrow",) + (("emphasis",) if row.emphasis else ())
                values = (
                    row.label,
                    f"{row.nominal:.2f}",
                    f"{row.favourable:.2f}",
                    f"{row.adverse:.2f}",
                    f"{row.mean:.2f}",
                    f"{row.variance:.3f}",
                    row.distribution,
                )
                tol_table.insert("", "end", values=values, tags=tags)
                export_rows.append(
                    (
                        "Tolerances (rain faded)",
                        row.label,
                        f"nominal {values[1]} / fav {values[2]} / adv {values[3]}",
                        f"mean {values[4]} / var {values[5]} / {values[6]}",
                    )
                )
            tol_table.pack(fill=tk.X)

        return export_rows

    def _render_pfd_section(self, parent, budget, extra_info, export_rows) -> None:
        """Power flux density check table, with notes on what is missing."""

        pfd_frame = ttk.Frame(parent)
        pfd_frame.pack(fill=tk.X)

        pfd_rows = static_report.pfd_rows(budget, extra_info)
        pfd_table = ttk.Treeview(pfd_frame, columns=("param", "value"), show="headings", height=len(pfd_rows))
        pfd_table.heading("param", text="Parameter")
        pfd_table.heading("value", text="Value")
        pfd_table.column("param", anchor="w", width=340)
        pfd_table.column("value", anchor="center", width=140)
        _configure_result_tags(pfd_table)
        row_help: dict[str, str] = {}
        for i, row in enumerate(pfd_rows):
            status = row.status()
            text = _format_checked(row.value, status, digits=row.digits)
            tags = ["evenrow" if i % 2 == 0 else "oddrow"] + ([status] if status else [])
            iid = pfd_table.insert("", "end", values=(row.label, text), tags=tuple(tags))
            if row.label in _FIXED_LB_ROW_HELP:
                row_help[iid] = _FIXED_LB_ROW_HELP[row.label]
            export_rows.append(("Power Flux Density", row.label, text, ""))
        _TreeviewRowTooltip(pfd_table, row_help)
        pfd_table.pack(anchor="w")

        for note in static_report.pfd_notes(budget, extra_info):
            ttk.Label(pfd_frame, text=note, foreground=PALETTE["danger"], wraplength=720, justify="left").pack(
                anchor="w", pady=(6, 0)
            )

    def calculate_ul_link_budget(self):
        """Render an uplink-only link budget table for the selected contact window."""

        if self.analysis_needs_refresh:
            messagebox.showwarning("Warning", "Please refresh the analysis first.")
            return

        if self.uplink_needs_refresh:
            self.recalculate_link_budget()

        if self.current_table_df.empty:
            messagebox.showwarning("Warning", "Please select a contact window first.")
            return

        required_columns = {
            "Time (UTC)",
            "Elevation (°)",
            "Slant Range (km)",
            "UL Path Loss (dB)",
            "UL Atmospheric Att (dB)",
            "UL Rx Power (dBW)",
            "UL C/No (dBHz)",
            "UL Eb/No (dB)",
        }
        missing_columns = required_columns - set(self.current_table_df.columns)
        if missing_columns:
            messagebox.showerror(
                "Data Error",
                "Uplink link budget values are not available. Please run the analysis again.",
            )
            return

        if self.uplink_table_frame is None:
            messagebox.showerror("UI Error", "Uplink results area is not initialised.")
            return

        self.clear_uplink_table()

        display_columns = [
            "Time (UTC)",
            "Elevation (°)",
            "Slant Range (km)",
            "UL Path Loss (dB)",
            "UL Atmospheric Att (dB)",
            "UL Rx Power (dBW)",
            "UL C/No (dBHz)",
            "UL Eb/No (dB)",
        ]

        def _format_cell(value):
            return "" if pd.isna(value) else value

        table = ttk.Treeview(self.uplink_table_frame, columns=display_columns, show="headings")
        table.tag_configure("evenrow", background=PALETTE["row_even"])
        table.tag_configure("oddrow", background=PALETTE["row_odd"])
        for col in display_columns:
            table.heading(col, text=col)
            if col == "Time (UTC)":
                table.column(col, anchor="center", width=150)
            else:
                table.column(col, anchor="center", width=130)

        for i, (_, row) in enumerate(self.current_table_df.iterrows()):
            formatted_time = row["Time (UTC)"].strftime("%Y-%m-%d %H:%M:%S")
            values = [formatted_time] + [_format_cell(row.get(col, "")) for col in display_columns[1:]]
            table.insert("", "end", values=values, tags=("evenrow" if i % 2 == 0 else "oddrow",))

        table_vscroll = ttk.Scrollbar(self.uplink_table_frame, orient=tk.VERTICAL, command=table.yview)
        table_hscroll = ttk.Scrollbar(self.uplink_table_frame, orient=tk.HORIZONTAL, command=table.xview)
        table.configure(yscrollcommand=table_vscroll.set, xscrollcommand=table_hscroll.set)
        table_vscroll.pack(side=tk.RIGHT, fill=tk.Y)
        table_hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        table.pack(fill=tk.BOTH, expand=True)
        self.set_uplink_budget_fresh()

    # ------------------------------------------------------------------
    # Baseband helpers
    # ------------------------------------------------------------------

    def _update_baseband_info(
        self,
        bit_rate_entry: ttk.Entry,
        roll_off_entry: ttk.Entry,
        overhead_entry: ttk.Entry,
        spectral_eff_entry: ttk.Entry,
        info_var: tk.StringVar,
        channel_bw_display: tk.StringVar,
        label_prefix: str,
    ):
        """Compute and display baseband derived metrics for a given entry set."""

        prefix = f"{label_prefix} " if label_prefix else ""

        try:
            bit_rate_mbps = float(bit_rate_entry.get())
            roll_off = float(roll_off_entry.get())
            overhead = float(overhead_entry.get())
            spectral_eff = float(spectral_eff_entry.get())
            info_bit_rate_mbps = bit_rate_mbps / overhead if overhead != 0 else 0
            spectral_eff = spectral_eff if spectral_eff != 0 else 1
            channel_bw_mhz = bit_rate_mbps * (1 + roll_off) / spectral_eff
            info_var.set(f"{prefix}Info Bit Rate [Mbps]: {info_bit_rate_mbps:.3f}")
            channel_bw_display.set(f"{prefix}Channel BW [MHz]: {channel_bw_mhz:.3f}")
        except ValueError:
            info_var.set(f"{prefix}Info Bit Rate [Mbps]: N/A")
            channel_bw_display.set(f"{prefix}Channel BW [MHz]: N/A")

    def update_link_budget_derived(self, *args):
        self._update_baseband_info(
            self.bitrate_entry,
            self.rolloff_entry,
            self.overhead_entry,
            self.spectral_eff_entry,
            self.info_bitrate_var,
            self.channel_bw_var,
            "",
        )

        if self.info_bitrate_ul_var is not None and self.channel_bw_ul_var is not None and all(
            entry is not None
            for entry in (
                self.uplink_bitrate_entry,
                self.uplink_rolloff_entry,
                self.uplink_overhead_entry,
                self.uplink_spectral_eff_entry,
            )
        ):
            self._update_baseband_info(
                self.uplink_bitrate_entry,
                self.uplink_rolloff_entry,
                self.uplink_overhead_entry,
                self.uplink_spectral_eff_entry,
                self.info_bitrate_ul_var,
                self.channel_bw_ul_var,
                "UL",
            )

    def on_panel_entry_edited(self, entry: ttk.Entry) -> None:
        """React to an edit of a parameters panel field, wherever it was typed.

        Shared by the pass analysis tabs and the static link budget tabs
        (which edit the same fields), so both refresh the derived baseband values and
        flag the uplink budget for recalculation in the same way.
        """

        uplink_baseband = (
            self.uplink_bitrate_entry,
            self.uplink_rolloff_entry,
            self.uplink_overhead_entry,
            self.uplink_spectral_eff_entry,
        )
        uplink_rf = (
            self.eirp_gs_entry,
            self.uplink_freq_entry,
            self.gt_sat_entry,
            self.demod_loss_ul_entry,
            self.other_att_ul_entry,
        )
        downlink_baseband = (self.bitrate_entry, self.rolloff_entry, self.overhead_entry, self.spectral_eff_entry)
        if any(entry is e for e in uplink_baseband):
            self.on_uplink_baseband_change()
        elif any(entry is e for e in uplink_rf):
            self.set_uplink_budget_stale()
        elif any(entry is e for e in downlink_baseband):
            self.update_link_budget_derived()

    def on_uplink_baseband_change(self, event=None):
        """Refresh derived UL metrics and mark uplink budget for recalculation."""

        self.update_link_budget_derived()
        self.set_uplink_budget_stale()

    # ------------------------------------------------------------------
    # GUI Setup
    # ------------------------------------------------------------------

    def _build_ui(self):
        """Initialise the Tkinter widget tree."""
        global _SPLASH_ROOT, _SPLASH_PROGRESS
        if _SPLASH_ROOT is not None:
            # Reuse the splash's root window instead of opening a second one.
            self.root = _SPLASH_ROOT
            _SPLASH_PROGRESS.stop()
            for widget in self.root.winfo_children():
                widget.destroy()
            self.root.overrideredirect(False)
            self.root.attributes("-topmost", False)
            _SPLASH_ROOT = None
            _SPLASH_PROGRESS = None
        else:
            self.root = tk.Tk()
        self.root.title(f"Satellite Link Budget Tool {APP_VERSION}")
        self.root.geometry("1200x950")
        self.root.minsize(1000, 700)
        self.root.configure(bg=PALETTE["bg"])

        base_font_family = _pick_font(["Segoe UI", "Helvetica Neue", "Helvetica", "Arial"])
        default_font = (base_font_family, 10)
        bold_font = (base_font_family, 10, "bold")
        header_font = (base_font_family, 15, "bold")
        subheader_font = (base_font_family, 10)

        style = ttk.Style()
        style.theme_use("clam")

        style.configure(".", background=PALETTE["bg"], foreground=PALETTE["text"], font=default_font)
        style.configure("TFrame", background=PALETTE["bg"])
        style.configure("Header.TFrame", background=PALETTE["header_bg"])
        style.configure("TLabel", background=PALETTE["bg"], foreground=PALETTE["text"])
        style.configure(
            "Header.TLabel", background=PALETTE["header_bg"], foreground=PALETTE["header_fg"], font=header_font
        )
        style.configure(
            "HeaderSub.TLabel",
            background=PALETTE["header_bg"],
            foreground=PALETTE["header_fg_muted"],
            font=subheader_font,
        )

        # Panels as MATLAB uipanels: grey like the window, a thin border and a
        # dark bold title.
        style.configure(
            "TLabelframe",
            background=PALETTE["bg"],
            bordercolor=PALETTE["border"],
            relief="solid",
            borderwidth=1,
        )
        style.configure(
            "TLabelframe.Label",
            background=PALETTE["bg"],
            foreground=PALETTE["text"],
            font=bold_font,
        )

        style.configure(
            "TButton",
            background="#fafafa",
            foreground=PALETTE["text"],
            padding=(10, 5),
            relief="solid",
            borderwidth=1,
            bordercolor="#b4b4b4",
            lightcolor="#fafafa",
            darkcolor="#fafafa",
            font=default_font,
        )
        style.map(
            "TButton",
            background=[("active", PALETTE["accent_light"]), ("disabled", PALETTE["bg"])],
            foreground=[("disabled", PALETTE["text_muted"])],
            bordercolor=[("active", PALETTE["accent"]), ("!active", "#b4b4b4")],
            lightcolor=[("active", PALETTE["accent_light"])],
            darkcolor=[("active", PALETTE["accent_light"])],
        )

        # Button on the dark header (Exit).
        style.configure(
            "Header.TButton",
            background=PALETTE["header_hover"],
            foreground=PALETTE["header_fg"],
            bordercolor=PALETTE["header_fg_muted"],
            lightcolor=PALETTE["header_hover"],
            darkcolor=PALETTE["header_hover"],
            padding=(14, 5),
            font=bold_font,
        )
        style.map(
            "Header.TButton",
            background=[("active", PALETTE["danger"])],
            bordercolor=[("active", PALETTE["danger"])],
            lightcolor=[("active", PALETTE["danger"])],
            darkcolor=[("active", PALETTE["danger"])],
        )

        # The action that needs running (no results yet, or inputs changed).
        style.configure(
            "Accent.TButton",
            background=PALETTE["accent"],
            foreground="white",
            font=bold_font,
            padding=(10, 5),
            relief="flat",
            borderwidth=1,
            bordercolor=PALETTE["accent_dark"],
            lightcolor=PALETTE["accent"],
            darkcolor=PALETTE["accent"],
        )
        style.map(
            "Accent.TButton",
            background=[("active", PALETTE["accent_dark"]), ("disabled", PALETTE["bg"])],
            foreground=[("disabled", PALETTE["text_muted"])],
            lightcolor=[("active", PALETTE["accent_dark"])],
            darkcolor=[("active", PALETTE["accent_dark"])],
        )
        # Large toolstrip buttons (static link budget actions); the
        # "Tool.Accent" style inherits the Accent colours.
        style.configure("Tool.TButton", padding=(14, 8))
        style.configure("Tool.Accent.TButton", padding=(14, 8))

        style.configure("TNotebook", background=PALETTE["bg"], borderwidth=0)
        # Main window tabs, as the MATLAB toolstrip tabs: a light grey strip,
        # the selected tab opening onto the page below with blue text.
        style.configure(
            "Main.TNotebook",
            background=PALETTE["toolstrip"],
            borderwidth=0,
            tabmargins=(8, 5, 8, 0),
            bordercolor=PALETTE["border"],
            lightcolor=PALETTE["toolstrip"],
            darkcolor=PALETTE["toolstrip"],
        )
        style.configure(
            "Main.TNotebook.Tab",
            background=PALETTE["toolstrip"],
            foreground="#3c3c3c",
            padding=(18, 7),
            font=bold_font,
            borderwidth=1,
            bordercolor=PALETTE["toolstrip"],
            lightcolor=PALETTE["toolstrip"],
        )
        style.map(
            "Main.TNotebook.Tab",
            background=[("selected", PALETTE["bg"]), ("active", PALETTE["toolstrip_hover"])],
            foreground=[("selected", PALETTE["accent"]), ("active", PALETTE["text"])],
            lightcolor=[("selected", PALETTE["bg"])],
            bordercolor=[("selected", PALETTE["border"])],
        )
        # Notebook without tabs, for windows that switch pages from a sidebar.
        style.layout("Pages.TNotebook.Tab", [])
        style.configure(
            "Pages.TNotebook",
            background=PALETTE["bg"],
            borderwidth=0,
            tabmargins=0,
            padding=0,
            bordercolor=PALETTE["bg"],
            lightcolor=PALETTE["bg"],
            darkcolor=PALETTE["bg"],
        )
        style.configure("Toolbar.TFrame", background=PALETTE["surface"])
        style.configure("Toolbar.TLabel", background=PALETTE["surface"])
        style.configure("Card.TLabel", background=PALETTE["surface"])
        style.configure(
            "ToolCaption.TLabel",
            background=PALETTE["surface"],
            foreground=PALETTE["text_muted"],
            font=(base_font_family, 8),
        )
        style.configure(
            "PageTitle.TLabel", foreground=PALETTE["accent_dark"], font=(base_font_family, 15, "bold")
        )
        style.configure("PageSub.TLabel", foreground=PALETTE["text_muted"], font=(base_font_family, 10))
        style.configure("Section.TLabel", foreground=PALETTE["accent_dark"], font=(base_font_family, 11, "bold"))
        style.configure(
            "TNotebook.Tab",
            background=PALETTE["surface_alt"],
            foreground=PALETTE["text_muted"],
            padding=(18, 9),
            font=bold_font,
            borderwidth=0,
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", PALETTE["surface"])],
            foreground=[("selected", PALETTE["accent_dark"])],
        )

        style.configure(
            "TEntry", fieldbackground=PALETTE["surface"], bordercolor=PALETTE["border"], padding=4
        )
        style.map("TEntry", bordercolor=[("focus", PALETTE["accent"]), ("!focus", PALETTE["border"])])

        style.configure(
            "TCombobox", fieldbackground=PALETTE["surface"], bordercolor=PALETTE["border"], padding=4
        )
        style.map(
            "TCombobox",
            bordercolor=[("focus", PALETTE["accent"]), ("!focus", PALETTE["border"])],
            fieldbackground=[("readonly", PALETTE["surface"])],
        )

        style.configure(
            "Treeview",
            background=PALETTE["surface"],
            fieldbackground=PALETTE["surface"],
            foreground=PALETTE["text"],
            rowheight=24,
            bordercolor=PALETTE["border"],
            borderwidth=1,
            font=default_font,
        )
        # Light column headings, as in MATLAB tables.
        style.configure(
            "Treeview.Heading",
            background=PALETTE["toolstrip"],
            foreground=PALETTE["text"],
            font=bold_font,
            relief="flat",
            bordercolor=PALETTE["border"],
            lightcolor=PALETTE["toolstrip"],
            darkcolor=PALETTE["border"],
            padding=(6, 5),
        )
        style.map("Treeview.Heading", background=[("active", PALETTE["toolstrip_hover"])])
        style.map(
            "Treeview",
            background=[("selected", PALETTE["selection"])],
            foreground=[("selected", PALETTE["text"])],
        )

        style.configure(
            "TScrollbar",
            background=PALETTE["surface_alt"],
            troughcolor=PALETTE["bg"],
            bordercolor=PALETTE["bg"],
            arrowsize=12,
        )

        # --- Header bar ---
        header = ttk.Frame(self.root, style="Header.TFrame", padding=(20, 14))
        header.grid(row=0, column=0, sticky="ew")
        _header_logo(header, base_font_family).pack(side=tk.LEFT)
        title_box = ttk.Frame(header, style="Header.TFrame")
        title_box.pack(side=tk.LEFT, padx=(16, 0))
        ttk.Label(title_box, text="Satellite Link Budget Tool", style="Header.TLabel").pack(anchor="w")
        ttk.Label(
            title_box,
            text=f"Version {APP_VERSION}  ·  Orbit prediction & downlink/uplink budget analysis",
            style="HeaderSub.TLabel",
        ).pack(anchor="w")
        # Far right of the header: Exit (also Ctrl+Q and the window close button).
        exit_button = ttk.Button(header, text="Exit", command=self.exit_app, style="Header.TButton")
        exit_button.pack(side=tk.RIGHT, padx=(18, 0))
        _add_tooltip(exit_button, "Close the program (Ctrl+Q).")
        about_button = ttk.Button(header, text="About", command=self.show_about, style="Header.TButton")
        about_button.pack(side=tk.RIGHT, padx=(18, 0))
        _add_tooltip(about_button, "Version, author and licence of the program (F1).")
        self.root.bind_all("<F1>", lambda _e: self.show_about())
        self.root.protocol("WM_DELETE_WINDOW", self.exit_app)
        self.root.bind_all("<Control-q>", lambda _e: self.exit_app())
        # Right of the header: the licence. The project name is in the window title.
        header_right = ttk.Frame(header, style="Header.TFrame")
        header_right.pack(side=tk.RIGHT, fill=tk.Y)
        self.project_var = tk.StringVar(value="")
        self.norad_var = tk.StringVar(value="")

        def _show_project(*_args):
            project = self.project_var.get().strip()
            title = f"Satellite Link Budget Tool {APP_VERSION}"
            self.root.title(f"{title} – {project}" if project else title)

        self.project_var.trace_add("write", _show_project)
        license_text = self.license_status.summary()
        if license_text:
            ttk.Label(
                header_right,
                text=license_text,
                style="HeaderSub.TLabel",
                foreground=PALETTE["warning_on_dark"] if self.license_status.expiring_soon else None,
            ).pack(side=tk.TOP, anchor="e", expand=True)

        self.root.grid_rowconfigure(0, weight=0)
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # Main tabs: the orbit and station set-up, the TLE-based pass analyses
        # and the static (fixed-elevation) link budgets.
        self.main_tabs = ttk.Notebook(self.root, style="Main.TNotebook")
        self.main_tabs.grid(row=1, column=0, sticky="nsew")

        def _scroll_tab(text: str) -> ttk.Frame:
            """Add a vertically scrollable tab; return the frame for its content."""

            container, body, tab_canvas = self._make_scrollable_frame(self.main_tabs)
            self.main_tabs.add(container, text=text)
            self._wheel_targets[str(container)] = lambda: tab_canvas
            content = ttk.Frame(body, padding=15)
            content.pack(fill=tk.BOTH, expand=True)
            return content

        def _wheel_target():
            selected = self.main_tabs.select()
            target = self._wheel_targets.get(selected)
            return target() if target else None

        # Base of the mouse-wheel stack: never popped, so any scrollable
        # popup's wheel binding reverts to this one when the popup closes.
        self._bind_scroll_wheel(self.main_tabs, target=_wheel_target)

        # The pass selected in the Orbit & Passes tab, shown above the results
        # of the Downlink and Uplink tabs.
        self.selected_pass_var = tk.StringVar(value=NO_PASS_SELECTED)

        def _static_losses_frame(parent, direction, include_var, text_var):
            frame = ttk.LabelFrame(parent, text="Link Losses from the Static Link Budget", padding=10)
            frame.pack(fill=tk.X, pady=(0, 10), padx=5)
            check = ttk.Checkbutton(
                frame,
                text=f"Include the losses entered in the Static LB · {direction} tab",
                variable=include_var,
                command=self._on_pass_losses_changed,
            )
            check.pack(anchor="w")
            _add_tooltip(
                check,
                "The pass analysis then subtracts the same elevation-independent losses as the "
                "static link budget: ionospheric, polarisation mismatch, multipath and antenna "
                "pointing, adds the modulation degradation to the demodulator loss, and refers "
                "Eb/No to the transfer frame bit rate (formatting factor). Atmospheric "
                "attenuation, EIRP, G/T and the other values of this tab are its own.",
            )
            ttk.Label(
                frame, textvariable=text_var, foreground=PALETTE["text_muted"], wraplength=1000, justify="left"
            ).pack(anchor="w", pady=(4, 0))

        def _selected_pass_label(parent, padx=0):
            row = ttk.Frame(parent)
            row.pack(fill=tk.X, padx=padx, pady=(0, 8))
            ttk.Label(row, text="Pass:", font=bold_font).pack(side=tk.LEFT)
            ttk.Label(
                row, textvariable=self.selected_pass_var, foreground=PALETTE["accent_dark"], font=bold_font
            ).pack(side=tk.LEFT, padx=(6, 0))
            ttk.Label(
                row, text="(select it in the Orbit & Passes tab)", foreground=PALETTE["text_muted"]
            ).pack(side=tk.LEFT, padx=(8, 0))

        setup_tab = _scroll_tab("Orbit & Passes")
        downlink_tab = _scroll_tab("Downlink")
        uplink_tab = _scroll_tab("Uplink")
        main_frame = setup_tab

        # --- Project ---
        project_frame = ttk.LabelFrame(main_frame, text="Satellite & Project", padding=10)
        project_frame.pack(fill=tk.X, pady=5)
        project_frame.grid_columnconfigure(0, weight=1)
        self.satellite_frame_rows = {}

        # Satellite list (satellites.json); hidden when there is none.
        select_row = ttk.Frame(project_frame)
        select_row.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        combo_line = ttk.Frame(select_row)
        combo_line.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(combo_line, text="Satellite", width=12).pack(side=tk.LEFT, padx=5)
        self.satellite_combo = ttk.Combobox(combo_line, state="readonly", width=60)
        self.satellite_combo.pack(side=tk.LEFT, padx=(6, 12))
        self.satellite_combo.bind("<<ComboboxSelected>>", self._on_satellite_chosen)
        _add_tooltip(
            self.satellite_combo,
            "The satellites of satellites.json. With missions/<satellite>/parameters.json the link budget "
            "is available; without it, the contact analysis only.",
        )
        self.mission_status_var = tk.StringVar(value="")
        ttk.Label(
            select_row, textvariable=self.mission_status_var, foreground=PALETTE["text_muted"], wraplength=1100
        ).pack(side=tk.TOP, anchor="w", padx=5, pady=(4, 0))
        select_row.grid_remove()
        self.satellite_frame_rows["select"] = select_row

        # Custom orbit of a mission not launched yet.
        self.orbit_alt_var = tk.StringVar(value="628")
        self.orbit_ecc_var = tk.StringVar(value="0")
        self.orbit_argp_var = tk.StringVar(value="0")
        self.orbit_ltan_var = tk.StringVar(value="22:20")
        self.orbit_inc_var = tk.StringVar(value="")
        orbit_row = ttk.Frame(project_frame)
        orbit_row.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        ttk.Label(orbit_row, text="Custom orbit", width=12).pack(side=tk.LEFT, padx=5)
        for label, var, width, tip in (
            (
                "Altitude (km)",
                self.orbit_alt_var,
                7,
                "Mean semi-major axis minus the equatorial radius (6378.137 km).",
            ),
            ("Ecc.", self.orbit_ecc_var, 9, "Mean eccentricity (0 for a circular orbit)."),
            ("ω (°)", self.orbit_argp_var, 5, "Mean argument of perigee (90° for a frozen orbit)."),
            (
                "Inclination (°)",
                self.orbit_inc_var,
                7,
                "Empty: sun-synchronous orbit, inclination computed from the altitude.",
            ),
            (
                "LTAN (hh:mm)",
                self.orbit_ltan_var,
                6,
                "Local time of the ascending node. For a constellation, choose the orbital plane.",
            ),
        ):
            ttk.Label(orbit_row, text=label).pack(side=tk.LEFT, padx=(6, 4))
            if var is self.orbit_ltan_var:
                orbit_entry = ttk.Combobox(orbit_row, textvariable=var, width=width)
                orbit_entry.bind("<<ComboboxSelected>>", lambda _e: self.generate_custom_tle(show_errors=False))
                self.orbit_ltan_combo = orbit_entry
            else:
                orbit_entry = ttk.Entry(orbit_row, textvariable=var, width=width)
            orbit_entry.pack(side=tk.LEFT)
            orbit_entry.bind("<KeyRelease>", lambda _e: self.set_analysis_stale())
            _add_tooltip(orbit_entry, tip)
        ttk.Button(orbit_row, text="Generate TLE", command=self.generate_custom_tle).pack(side=tk.LEFT, padx=(12, 0))
        orbit_row.grid_remove()
        self.satellite_frame_rows["orbit"] = orbit_row

        name_row = ttk.Frame(project_frame)
        name_row.grid(row=2, column=0, sticky="ew")
        ttk.Label(name_row, text="Project Name", width=12).pack(side=tk.LEFT, padx=5)
        project_entry = ttk.Entry(name_row, textvariable=self.project_var, width=40)
        project_entry.pack(side=tk.LEFT, padx=(6, 12))
        ttk.Label(
            name_row,
            text="Shown in the window title and in the reports; saved in the parameters file.",
            foreground=PALETTE["text_muted"],
        ).pack(side=tk.LEFT)

        # --- TLE Section ---
        tle_frame = ttk.LabelFrame(main_frame, text="TLE Parameters", padding=10)
        tle_frame.pack(fill=tk.X, pady=5)
        ttk.Label(tle_frame, text="TLE Line 1").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.tle1_entry = ttk.Entry(tle_frame, width=80)
        self.tle1_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        self.tle1_entry.bind("<KeyRelease>", lambda event: self.set_analysis_stale())
        ttk.Label(tle_frame, text="TLE Line 2").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.tle2_entry = ttk.Entry(tle_frame, width=80)
        self.tle2_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        self.tle2_entry.bind("<KeyRelease>", lambda event: self.set_analysis_stale())
        tle_buttons = ttk.Frame(tle_frame)
        tle_buttons.grid(row=2, column=0, columnspan=2, sticky="w", padx=5, pady=(6, 0))
        ttk.Button(tle_buttons, text="Load TLE from file", command=self.load_tle_from_file).pack(side=tk.LEFT)
        ttk.Label(tle_buttons, text="NORAD ID").pack(side=tk.LEFT, padx=(24, 6))
        norad_entry = ttk.Entry(tle_buttons, textvariable=self.norad_var, width=8)
        norad_entry.pack(side=tk.LEFT)
        self.tle_download_button = ttk.Button(tle_buttons, text="Download TLE", command=self.download_tle)
        self.tle_download_button.pack(side=tk.LEFT, padx=(6, 0))
        _add_tooltip(
            self.tle_download_button,
            "Download the latest TLE of this NORAD ID from SatNOGS (or CelesTrak) and save it in tle.txt.",
        )
        self.tle_file_var = tk.StringVar(value="TLE: none loaded")
        self.tle_info_label = ttk.Label(
            tle_frame, textvariable=self.tle_file_var, foreground=PALETTE["text_muted"]
        )
        self.tle_info_label.grid(row=3, column=0, columnspan=2, sticky="w", padx=5, pady=(4, 0))
        tle_frame.grid_columnconfigure(1, weight=1)

        # --- Observation Settings Section ---
        obs_frame = ttk.LabelFrame(main_frame, text="Observation Settings", padding=10)
        obs_frame.pack(fill=tk.X, pady=5)
        obs_frame.grid_columnconfigure(0, weight=1)

        # Row 1: what pass to analyse (date + ground station).
        pass_row = ttk.Frame(obs_frame)
        pass_row.grid(row=0, column=0, sticky="ew", padx=5)

        ttk.Label(pass_row, text="Date (YYYY-MM-DD)").pack(side=tk.LEFT)
        self.date_entry = ttk.Entry(pass_row, width=15)
        self.date_entry.pack(side=tk.LEFT, padx=(6, 24))
        self.date_entry.insert(0, datetime.now(timezone.utc).strftime("%Y-%m-%d"))
        self.date_entry.bind("<KeyRelease>", lambda event: self.set_analysis_stale())

        ttk.Label(pass_row, text="Ground Station").pack(side=tk.LEFT)
        initial_stations = list(calculations.GROUND_STATIONS.keys())
        self.gs_var = tk.StringVar(value=initial_stations[0] if initial_stations else "")
        self.gs_menu = ttk.Combobox(
            pass_row, textvariable=self.gs_var, values=initial_stations, state="readonly", width=15
        )
        self.gs_menu.pack(side=tk.LEFT, padx=(6, 12))
        self.gs_menu.bind("<<ComboboxSelected>>", lambda event: self.set_analysis_stale())

        load_gs_btn = ttk.Button(
            pass_row, text="Load Ground Stations...", command=self.load_ground_stations_from_file
        )
        load_gs_btn.pack(side=tk.LEFT)
        _add_tooltip(
            load_gs_btn,
            "Load a custom ground station catalogue (name, latitude, longitude, altitude) from a "
            ".txt/.csv file. Replaces the choices in the Ground Station dropdown above.",
        )

        self.gs_file_var = tk.StringVar(value="Ground stations: none loaded")
        ttk.Label(obs_frame, textvariable=self.gs_file_var, foreground=PALETTE["text_muted"]).grid(
            row=1, column=0, sticky="w", padx=5, pady=(4, 10)
        )

        ttk.Separator(obs_frame, orient="horizontal").grid(row=2, column=0, sticky="ew", pady=(0, 10))

        # Row 2: bulk-import link budget parameters from a JSON file.
        params_row = ttk.Frame(obs_frame)
        params_row.grid(row=3, column=0, sticky="ew", padx=5)
        load_params_btn = ttk.Button(
            params_row,
            text="Import Link Budget Parameters (.json)...",
            command=self.load_parameters_from_file,
        )
        load_params_btn.pack(side=tk.LEFT)
        params_hint = (
            "Fills EIRP, G/T, bit rate, roll-off and the other downlink/uplink fields below from a "
            "saved JSON file, instead of typing them in one by one."
        )
        _add_tooltip(load_params_btn, f"{params_hint} See the README for the expected file format.")
        ttk.Label(params_row, text=params_hint, foreground=PALETTE["text_muted"], wraplength=620).pack(
            side=tk.LEFT, padx=(12, 0)
        )

        self.param_file_var = tk.StringVar(value="Parameters: none loaded")
        ttk.Label(obs_frame, textvariable=self.param_file_var, foreground=PALETTE["text_muted"]).grid(
            row=4, column=0, sticky="w", padx=5, pady=(4, 0)
        )

        # --- Pass analysis (shared by the Downlink and Uplink tabs) ---
        passes_frame = ttk.LabelFrame(main_frame, text="Pass Analysis", padding=10)
        passes_frame.pack(fill=tk.X, pady=5)
        passes_row = ttk.Frame(passes_frame)
        passes_row.pack(fill=tk.X, padx=5)
        self.start_refresh_button = ttk.Button(
            passes_row,
            text="Start/Refresh Analysis",
            command=self.run_analysis,
            style="Accent.TButton",
        )
        self.start_refresh_button.pack(side=tk.LEFT)
        _add_tooltip(
            self.start_refresh_button,
            "Re-propagates the orbit for the full day (TLE, date, ground station) and finds all "
            "contact windows. Turns blue whenever those inputs change and a fresh run is needed.",
        )
        ttk.Label(
            passes_row,
            text=(
                "Finds the passes over the ground station on the chosen date, with the downlink and "
                "uplink link budgets along each pass. Select a pass below: the Downlink and Uplink "
                "tabs show its results."
            ),
            foreground=PALETTE["text_muted"],
            wraplength=760,
            justify="left",
        ).pack(side=tk.LEFT, padx=(14, 0))

        antenna_row = ttk.Frame(passes_frame)
        antenna_row.pack(fill=tk.X, padx=5, pady=(10, 0))
        ttk.Label(antenna_row, text="Satellite antenna pattern (gain vs angle from nadir):").pack(side=tk.LEFT)
        _info_icon(
            antenna_row,
            "Satellite antenna gain versus the angle from nadir (boresight), applied at every time "
            "step in both directions. Downlink: EIRP = Tx chain of the Static LB · Downlink tab + "
            "the pattern gain towards the station (or EIRP SAT at the pattern peak when that tab "
            "uses an EIRP value). Uplink: G/T SAT at the pattern peak minus the roll-off. The "
            "ground station antenna tracks the satellite: no pointing loss. At startup the pattern "
            f"is read from '{ANTENNA_PATTERN_FILE}' (next to the program), if present.",
        ).pack(side=tk.LEFT, padx=(4, 0))
        load_antenna_btn = ttk.Button(
            antenna_row, text="Load Antenna Pattern...", command=self.load_antenna_pattern_file
        )
        load_antenna_btn.pack(side=tk.LEFT, padx=(12, 6))
        _add_tooltip(
            load_antenna_btn,
            "Load the satellite antenna gain pattern from a .csv/.txt file: two columns, angle "
            "from the antenna boresight (nadir) in degrees and gain in dBi.",
        )
        ttk.Button(antenna_row, text="Show Antenna Gain", command=self.show_antenna_pattern).pack(side=tk.LEFT)
        self.pattern_file_var = tk.StringVar(value="Pattern: built-in example (no file loaded)")
        ttk.Label(passes_frame, textvariable=self.pattern_file_var, foreground=PALETTE["text_muted"]).pack(
            anchor="w", padx=5, pady=(4, 0)
        )

        contact_frame = ttk.Frame(passes_frame)
        contact_frame.pack(fill=tk.X, padx=5, pady=(10, 0))
        ttk.Label(contact_frame, text="Contact Windows (UTC Time)", font=bold_font).pack(anchor="w", pady=(0, 4))
        self.contact_mode_var = tk.StringVar(value="")
        ttk.Label(contact_frame, textvariable=self.contact_mode_var, foreground=PALETTE["warning"]).pack(
            anchor="w", pady=(0, 4)
        )
        self.contact_listbox = tk.Listbox(
            contact_frame,
            height=6,
            exportselection=False,
            bg=PALETTE["surface"],
            fg=PALETTE["text"],
            selectbackground=PALETTE["accent"],
            selectforeground="white",
            activestyle="none",
            relief="flat",
            borderwidth=0,
            highlightthickness=1,
            highlightbackground=PALETTE["border"],
            highlightcolor=PALETTE["accent"],
            font=default_font,
        )
        self.contact_listbox.pack(fill=tk.X, expand=True)
        self.contact_listbox.bind("<<ListboxSelect>>", self.on_contact_select)


        param_container = ttk.Frame(downlink_tab)
        param_container.pack(fill=tk.X, expand=False, padx=5, pady=5)
        param_container.grid_columnconfigure(0, weight=1)
        param_container.grid_columnconfigure(1, weight=1)
        param_container.grid_columnconfigure(2, weight=1)
        param_container.grid_rowconfigure(0, weight=1)
        param_container.grid_rowconfigure(1, weight=1)

        # Satellite Parameters frame
        sat_frame = ttk.LabelFrame(param_container, text="Satellite Parameters", padding=10)
        sat_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 10), pady=5)
        ttk.Label(sat_frame, text="EIRP SAT [dBW]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.eirp_sat_entry = ttk.Entry(sat_frame, width=15)
        self.eirp_sat_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(sat_frame, text="Frequency [GHz]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.freq_entry = ttk.Entry(sat_frame, width=15)
        self.freq_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(sat_frame, text="C/Io [dBHz]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.CIo_entry = ttk.Entry(sat_frame, width=15)
        self.CIo_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
        _add_tooltip(
            self.CIo_entry,
            "Carrier to interference density ratio C/Io, combined with C/No into C/(No+Io). "
            "Leave empty when there is no interferer: the link budget then uses C/No alone.",
        )
        sat_frame.grid_columnconfigure(1, weight=1)

        # Ground Station Parameters frame
        gs_frame = ttk.LabelFrame(param_container, text="Ground Station Parameters", padding=10)
        gs_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 10), pady=5)
        ttk.Label(gs_frame, text="G/T GS [dB/K]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.gt_gs_entry = ttk.Entry(gs_frame, width=15)
        self.gt_gs_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(gs_frame, text="Antenna Diameter GS [m]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.d_gs_entry = ttk.Entry(gs_frame, width=15)
        self.d_gs_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        gs_frame.grid_columnconfigure(1, weight=1)

        # Atmospheric Attenuations frame
        atm_frame = ttk.LabelFrame(param_container, text="Atmospheric Attenuations", padding=10)
        atm_frame.grid(row=0, column=2, sticky="nsew", pady=5)
        ttk.Label(atm_frame, text="Link Availability [%]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.LA_entry = ttk.Entry(atm_frame, width=15)
        self.LA_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(atm_frame, text="Other Attenuations [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.other_att_entry = ttk.Entry(atm_frame, width=15)
        self.other_att_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        # Off by default, as in the static link budget: the ITU-R P.618 model
        # is given for elevations of at least 5° and about 4-20 GHz.
        self.scint_var = tk.BooleanVar(value=False)
        scint_check = ttk.Checkbutton(
            atm_frame,
            text="Tropospheric scintillation",
            variable=self.scint_var,
            command=self.set_analysis_stale,
        )
        scint_check.grid(row=2, column=0, columnspan=2, sticky="w", padx=5, pady=(8, 0))
        _add_tooltip(
            scint_check,
            "ITU-R P.618 tropospheric scintillation fading, applied to both downlink and uplink. "
            "Off by default: the model is given for elevations of at least 5° and about 4–20 GHz "
            "and matters mainly in Ku/Ka-band; below that (e.g. L/S band at low elevation) it is "
            "an extrapolation. Ionospheric scintillation (the relevant one at L/S band) comes from "
            "the Static LB tabs (S4 index), with the other link losses.",
        )
        self.atm_label_var = tk.StringVar(
            value=f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}° El: N/A"
        )
        ttk.Label(atm_frame, textvariable=self.atm_label_var).grid(row=3, column=0, columnspan=2, sticky="w", pady=(10, 0), padx=5)
        atm_frame.grid_columnconfigure(1, weight=1)

        # Baseband Parameters frame
        baseband_frame = ttk.LabelFrame(param_container, text="Baseband Parameters", padding=10)
        baseband_frame.grid(row=1, column=0, columnspan=3, sticky="nsew", pady=(0, 5))
        ttk.Label(baseband_frame, text="Coded Bit Rate [Mbps]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.bitrate_entry = ttk.Entry(baseband_frame, width=15)
        self.bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Roll-off Factor").grid(row=0, column=2, sticky="w", padx=5, pady=2)
        self.rolloff_entry = ttk.Entry(baseband_frame, width=15)
        self.rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Demodulator Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.demod_loss_entry = ttk.Entry(baseband_frame, width=15)
        self.demod_loss_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Overhead (coded / info)").grid(row=1, column=2, sticky="w", padx=5, pady=2)
        self.overhead_entry = ttk.Entry(baseband_frame, width=15)
        self.overhead_entry.grid(row=1, column=3, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Spectral Efficiency [b/s/Hz]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.spectral_eff_entry = ttk.Entry(baseband_frame, width=15)
        self.spectral_eff_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
        self.info_bitrate_var = tk.StringVar(value="Info Bit Rate [Mbps]: N/A")
        self.channel_bw_var = tk.StringVar(value="Channel BW [MHz]: N/A")
        ttk.Label(baseband_frame, textvariable=self.info_bitrate_var).grid(row=3, column=0, columnspan=2, sticky="w", pady=(5, 0), padx=5)
        ttk.Label(baseband_frame, textvariable=self.channel_bw_var).grid(row=3, column=2, columnspan=2, sticky="w", padx=5)
        baseband_frame.grid_columnconfigure(1, weight=1)
        baseband_frame.grid_columnconfigure(3, weight=1)
        for entry in (self.bitrate_entry, self.rolloff_entry, self.overhead_entry, self.spectral_eff_entry):
            entry.bind("<KeyRelease>", lambda event, edited=entry: self.on_panel_entry_edited(edited))
        self.update_link_budget_derived()

        self.dl_losses_var = tk.StringVar(value="")
        self.include_static_dl_var = tk.BooleanVar(value=True)
        _static_losses_frame(downlink_tab, "Downlink", self.include_static_dl_var, self.dl_losses_var)

        downlink_btn_frame = ttk.LabelFrame(
            downlink_tab, text="Downlink Actions", padding=10
        )
        downlink_btn_frame.pack(fill=tk.X, pady=(0, 10), padx=5)
        self.recalc_button = ttk.Button(
            downlink_btn_frame,
            text="Calculate DL Link Budget",
            command=self.recalculate_link_budget,
            state="disabled",
        )
        self.recalc_button.pack(side=tk.LEFT, padx=5)
        _add_tooltip(
            self.recalc_button,
            "Re-computes only the downlink link budget (EIRP, G/T, attenuations, bit rate, ...) for "
            "the passes already found in the Orbit & Passes tab, without re-running the orbit "
            "propagation.",
        )

        downlink_results_container = ttk.Frame(downlink_tab)
        downlink_results_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))

        ttk.Button(downlink_btn_frame, text="Export Table CSV", command=self.export_table_csv).pack(
            side=tk.RIGHT, padx=5
        )

        _selected_pass_label(downlink_results_container)

        # Fixed height: the plots are often drawn while this tab is hidden
        # (a pass selected in the Orbit & Passes tab), and a canvas sized
        # by its content would collapse there.
        self.plot_frame = ttk.Frame(downlink_results_container, height=400)
        self.plot_frame.pack(fill=tk.X, pady=(0, 5))
        self.plot_frame.pack_propagate(False)

        downlink_results_frame = ttk.LabelFrame(
            downlink_results_container, text="Downlink Results", padding=10
        )
        downlink_results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))
        self.table_frame = ttk.Frame(downlink_results_frame)
        self.table_frame.pack(fill=tk.BOTH, expand=True)

        # Uplink Parameters frame
        uplink_container = ttk.Frame(uplink_tab)
        uplink_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        uplink_container.grid_columnconfigure(0, weight=1)
        uplink_container.grid_columnconfigure(1, weight=1)
        uplink_container.grid_columnconfigure(2, weight=1)
        uplink_container.grid_rowconfigure(0, weight=1)
        uplink_container.grid_rowconfigure(1, weight=1)

        uplink_gs_frame = ttk.LabelFrame(uplink_container, text="Ground Station Parameters", padding=10)
        uplink_gs_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 10), pady=5)
        ttk.Label(uplink_gs_frame, text="EIRP GS [dBW]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.eirp_gs_entry = ttk.Entry(uplink_gs_frame, width=15)
        self.eirp_gs_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_gs_frame, text="Uplink Frequency [GHz]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.uplink_freq_entry = ttk.Entry(uplink_gs_frame, width=15)
        self.uplink_freq_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        uplink_gs_frame.grid_columnconfigure(1, weight=1)

        uplink_sat_frame = ttk.LabelFrame(uplink_container, text="Satellite Parameters", padding=10)
        uplink_sat_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 10), pady=5)
        ttk.Label(uplink_sat_frame, text="G/T SAT [dB/K]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.gt_sat_entry = ttk.Entry(uplink_sat_frame, width=15)
        self.gt_sat_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_sat_frame, text="Sat Demod Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.demod_loss_ul_entry = ttk.Entry(uplink_sat_frame, width=15)
        self.demod_loss_ul_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        uplink_sat_frame.grid_columnconfigure(1, weight=1)

        uplink_atm_frame = ttk.LabelFrame(uplink_container, text="Uplink Attenuations", padding=10)
        uplink_atm_frame.grid(row=0, column=2, sticky="nsew", pady=5)
        ttk.Label(uplink_atm_frame, text="Uplink Other Attenuations [dB]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.other_att_ul_entry = ttk.Entry(uplink_atm_frame, width=15)
        self.other_att_ul_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        uplink_atm_frame.grid_columnconfigure(1, weight=1)

        uplink_baseband_frame = ttk.LabelFrame(
            uplink_container, text="Uplink Baseband Parameters", padding=10
        )
        uplink_baseband_frame.grid(row=1, column=0, columnspan=3, sticky="nsew", pady=(0, 5))
        ttk.Label(uplink_baseband_frame, text="UL Coded Bit Rate [Mbps]").grid(
            row=0, column=0, sticky="w", padx=5, pady=2
        )
        self.uplink_bitrate_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_baseband_frame, text="UL Roll-off Factor").grid(
            row=0, column=2, sticky="w", padx=5, pady=2
        )
        self.uplink_rolloff_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_baseband_frame, text="UL Overhead (coded / info)").grid(
            row=1, column=0, sticky="w", padx=5, pady=2
        )
        self.uplink_overhead_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_overhead_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_baseband_frame, text="UL Spectral Efficiency [b/s/Hz]").grid(
            row=1, column=2, sticky="w", padx=5, pady=2
        )
        self.uplink_spectral_eff_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_spectral_eff_entry.grid(row=1, column=3, sticky="ew", padx=5, pady=2)
        self.info_bitrate_ul_var = tk.StringVar(value="UL Info Bit Rate [Mbps]: N/A")
        self.channel_bw_ul_var = tk.StringVar(value="UL Channel BW [MHz]: N/A")
        ttk.Label(uplink_baseband_frame, textvariable=self.info_bitrate_ul_var).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(5, 0), padx=5
        )
        ttk.Label(uplink_baseband_frame, textvariable=self.channel_bw_ul_var).grid(
            row=2, column=2, columnspan=2, sticky="w", padx=5
        )

        uplink_baseband_frame.grid_columnconfigure(1, weight=1)
        uplink_baseband_frame.grid_columnconfigure(3, weight=1)

        for entry in [
            self.eirp_gs_entry,
            self.uplink_freq_entry,
            self.gt_sat_entry,
            self.demod_loss_ul_entry,
            self.other_att_ul_entry,
        ]:
            entry.bind("<KeyRelease>", lambda event, edited=entry: self.on_panel_entry_edited(edited))

        for entry in [
            self.uplink_bitrate_entry,
            self.uplink_rolloff_entry,
            self.uplink_overhead_entry,
            self.uplink_spectral_eff_entry,
        ]:
            entry.bind("<KeyRelease>", lambda event, edited=entry: self.on_panel_entry_edited(edited))

        self.ul_losses_var = tk.StringVar(value="")
        self.include_static_ul_var = tk.BooleanVar(value=True)
        _static_losses_frame(uplink_tab, "Uplink", self.include_static_ul_var, self.ul_losses_var)

        uplink_actions_frame = ttk.LabelFrame(
            uplink_tab, text="Uplink Actions", padding=10
        )
        uplink_actions_frame.pack(fill=tk.X, pady=10, padx=5)
        self.uplink_recalc_button = ttk.Button(
            uplink_actions_frame,
            text="Calculate UL Link Budget",
            command=self.calculate_ul_link_budget,
        )
        self.uplink_recalc_button.pack(side=tk.LEFT, padx=5)

        _selected_pass_label(uplink_tab, padx=5)

        uplink_results_frame = ttk.LabelFrame(uplink_tab, text="Uplink Results", padding=10)
        uplink_results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        self.uplink_table_frame = ttk.Frame(uplink_results_frame)
        self.uplink_table_frame.pack(fill=tk.BOTH, expand=True)

        # Static (fixed-elevation) link budgets, after the panels whose
        # fields they share.
        self.static_dl_tab = ttk.Frame(self.main_tabs)
        self.main_tabs.add(self.static_dl_tab, text="Static LB · Downlink")
        self._build_static_downlink_view(self.static_dl_tab)
        self.static_ul_tab = ttk.Frame(self.main_tabs)
        self.main_tabs.add(self.static_ul_tab, text="Static LB · Uplink")
        self._build_static_uplink_view(self.static_ul_tab)

        # Window icon, also inherited by every popup window.
        icon_path = _asset_path(ICON_FILE)
        if icon_path and sys.platform == "win32":
            try:
                self.root.iconbitmap(default=icon_path)
            except tk.TclError:
                pass
        icon_image = _load_image(ICON_PNG_FILE, master=self.root)
        if icon_image is not None:
            self._icon_image = icon_image
            self.root.iconphoto(True, icon_image)

        self._load_startup_defaults()
        # Show the losses the pass analysis takes from the static tabs, kept
        # up to date whenever a tab is opened.
        self._pass_loss_settings()
        self.main_tabs.bind("<<NotebookTabChanged>>", lambda _e: self._pass_loss_settings(), add="+")


def setup_gui():
    """Initialise and launch the Tkinter graphical interface.

    Kept as a module-level, zero-argument entry point for backward
    compatibility with ``__main__.py`` and ``__init__.py``.
    """
    app = LinkBudgetApp()
    app.run()
    return app


if __name__ == "__main__":
    setup_gui()
