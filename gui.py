"""Tkinter GUI for the satellite link budget tool."""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import tkinter.font as tkfont

# ---------------------------------------------------------------------------
# Visual theme
# ---------------------------------------------------------------------------
# Defined up here, before the slower imports below, because the startup
# splash screen needs it immediately.
PALETTE = {
    "bg": "#eef1f6",
    "surface": "#ffffff",
    "surface_alt": "#f7f9fc",
    "border": "#d7dde6",
    "header_bg": "#0f2440",
    "header_fg": "#f5f8ff",
    "header_fg_muted": "#a9b8d1",
    "accent": "#2f6fed",
    "accent_dark": "#2354bd",
    "accent_light": "#e8f0fe",
    "text": "#1f2733",
    "text_muted": "#5b6472",
    "danger": "#d64545",
    "danger_dark": "#b23434",
    "row_odd": "#ffffff",
    "row_even": "#f3f6fb",
    "selection": "#cfe0fd",
}


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

    width, height = 420, 170
    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()
    root.geometry(f"{width}x{height}+{(screen_w - width) // 2}+{(screen_h - height) // 2}")

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
        text="Loading, please wait...",
        bg=PALETTE["header_bg"],
        fg=PALETTE["header_fg_muted"],
        font=(font_family, 9),
    ).pack(pady=(4, 14))

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

from datetime import datetime, timedelta, timezone

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MinuteLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import pandas as pd

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
from skyfield.api import load, EarthSatellite, wgs84

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import astropy.units as u
import os
import sys
import json

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)
import calculations

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
    atmospheric_attenuation,
    prepare_topocentric_data,
    load_antenna_pattern,
    reload_ground_stations,
    MIN_ELEVATION_DEG,
)

_splash_tick(_SPLASH_ROOT, _SPLASH_PROGRESS)


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
        self.param_file_var: tk.StringVar | None = None
        self.uplink_bitrate_entry: ttk.Entry | None = None
        self.uplink_rolloff_entry: ttk.Entry | None = None
        self.uplink_overhead_entry: ttk.Entry | None = None
        self.uplink_spectral_eff_entry: ttk.Entry | None = None
        self.info_bitrate_ul_var: tk.StringVar | None = None
        self.channel_bw_ul_var: tk.StringVar | None = None
        self.uplink_table_frame: ttk.Frame | None = None
        self.uplink_recalc_button: ttk.Button | None = None

        self._build_ui()

    def run(self):
        """Enter the Tk event loop; blocks until the window is closed."""
        self.root.mainloop()

    # ------------------------------------------------------------------
    # File loading callbacks
    # ------------------------------------------------------------------

    def load_tle_from_file(self):
        """Load TLE lines from a text file and populate the GUI fields."""

        file_path = filedialog.askopenfilename(
            title="Select TLE file",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
        )
        if not file_path:
            return

        try:
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
        except Exception as exc:
            messagebox.showerror("TLE Error", f"Failed to load TLE: {exc}")
            return

        self.tle1_entry.delete(0, tk.END)
        self.tle1_entry.insert(0, tle1)
        self.tle2_entry.delete(0, tk.END)
        self.tle2_entry.insert(0, tle2)
        self.set_analysis_stale()

    def load_parameters_from_file(self):
        """Load all configurable parameters from a JSON file.

        The JSON structure should provide keys such as ``eirp_sat_dbw`` or
        ``frequency_ghz``. See the README for the full schema and example.
        Missing keys are ignored so the user can provide only the fields
        they need.
        """

        file_path = filedialog.askopenfilename(
            title="Select Parameters File",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
        )
        if not file_path:
            return

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception as exc:
            messagebox.showerror("Parameters", f"Unable to load parameters: {exc}")
            return

        field_map = {
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

        for key, entry in field_map.items():
            if key in payload and entry is not None:
                entry.delete(0, tk.END)
                entry.insert(0, str(payload[key]))

        if self.param_file_var is not None:
            self.param_file_var.set(f"Parameters: {os.path.abspath(file_path)}")

        self.update_link_budget_derived()
        self.set_analysis_stale()

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
            self.start_refresh_button.config(style="Red.TButton")
            self.recalc_button.config(state="disabled")
            self.clear_plot_and_table()
            self.clear_uplink_table()
            self.contact_listbox.delete(0, tk.END)

    def set_uplink_budget_stale(self):
        """Mark the uplink-only link budget as stale."""

        if not self.uplink_needs_refresh:
            self.uplink_needs_refresh = True
            self.clear_uplink_table()
            if self.uplink_recalc_button is not None:
                self.uplink_recalc_button.config(style="Red.TButton")

    def set_uplink_budget_fresh(self):
        """Reset the uplink-only link budget stale flag."""

        if self.uplink_needs_refresh:
            self.uplink_needs_refresh = False
            if self.uplink_recalc_button is not None:
                self.uplink_recalc_button.config(style="TButton")

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
            )
            rounding_rules = {"Atmospheric Att (dB)": 3, "UL Atmospheric Att (dB)": 3}
            for key in [
                "Elevation (°)",
                "Slant Range (km)",
                "Path Loss (dB)",
                "Atmospheric Att (dB)",
                "Pointing Loss (dB)",
                "Off Boresight Angle (°)",
                "Rx Power (dBW)",
                "C/(No+Io) (dBHz)",
                "Eb/No (dB)",
                "Doppler Shift (kHz)",
                "UL Path Loss (dB)",
                "UL Atmospheric Att (dB)",
                "UL Pointing Loss (dB)",
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

        self.analysis_needs_refresh = False
        self.start_refresh_button.config(style="TButton")
        self.recalc_button.config(state="!disabled")

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
            )

            rounding_rules = {"Atmospheric Att (dB)": 3, "UL Atmospheric Att (dB)": 3}
            for key in [
                "Elevation (°)",
                "Slant Range (km)",
                "Path Loss (dB)",
                "Atmospheric Att (dB)",
                "Pointing Loss (dB)",
                "Off Boresight Angle (°)",
                "Rx Power (dBW)",
                "C/(No+Io) (dBHz)",
                "Eb/No (dB)",
                "Doppler Shift (kHz)",
                "UL Path Loss (dB)",
                "UL Atmospheric Att (dB)",
                "UL Pointing Loss (dB)",
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
            self.contact_listbox.insert(
                tk.END,
                f"Contact {i+1}: {start.strftime('%H:%M:%S')} - {end.strftime('%H:%M:%S')} ({pass_direction})",
            )

    def on_contact_select(self, event):
        selection = self.contact_listbox.curselection()
        if not selection:
            self.clear_plot_and_table()
            self.clear_uplink_table()
            return

        idx = selection[0]
        start, end, _ = self.contact_windows[idx]
        mask = (self.df_all["Time (UTC)"] >= start) & (self.df_all["Time (UTC)"] <= end)
        df_pass = self.df_all[mask].copy()
        self.current_table_df = df_pass

        self.clear_uplink_table()
        self.clear_plot_and_table()

        fig, (ax1, ax2, ax3) = plt.subplots(
            1,
            3,
            figsize=(15, 4),
            gridspec_kw={"width_ratios": [1, 1, 1]},
        )

        tk_bg_color = self.root.cget("bg")
        try:
            r, g, b = self.root.winfo_rgb(tk_bg_color)
            mpl_bg_color = f"#{r // 256:02x}{g // 256:02x}{b // 256:02x}"
        except tk.TclError:
            mpl_bg_color = "#F0F0F0"
        fig.patch.set_facecolor(mpl_bg_color)

        ax1.plot(df_pass["Time (UTC)"], df_pass["Elevation (°)"], label="Elevation (°)", color="tab:blue")
        ax1.set_ylabel("Elevation (°)", color="tab:blue", fontsize=8)
        ax1.tick_params(axis="both", labelcolor="tab:blue", labelsize=8)
        ax1.set_xlabel("Time (UTC)", fontsize=8)
        ax1.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.grid(True, linestyle=":", alpha=0.7)
        ax1.set_title("Elevation Angle")
        fig.autofmt_xdate(rotation=30)

        ax2_1 = ax2
        ax2_2 = ax2.twinx()
        ax2_1.plot(df_pass["Time (UTC)"], df_pass["Eb/No (dB)"], label="Eb/No (dB)", color="tab:red")
        ax2_1.set_ylabel("Eb/No (dB)", color="tab:red", fontsize=8)
        ax2_1.tick_params(axis="both", labelcolor="tab:red", labelsize=8)
        ax2_2.plot(
            df_pass["Time (UTC)"],
            df_pass["C/(No+Io) (dBHz)"],
            label="C/(No+Io) (dBHz)",
            color="tab:green",
            linestyle="--",
        )
        ax2_2.set_ylabel("C/(No+Io) (dBHz)", color="tab:green", fontsize=8)
        ax2_2.tick_params(axis="both", labelcolor="tab:green", labelsize=8)
        ax2.set_xlabel("Time (UTC)", fontsize=8)
        ax2.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.grid(True, linestyle=":", alpha=0.7)
        ax2.set_title("Link Quality (Eb/No & C/(No+Io))")
        lines, labels = ax2_1.get_legend_handles_labels()
        lines2, labels2 = ax2_2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc="best")

        ax3.plot(df_pass["Time (UTC)"], df_pass["Doppler Shift (kHz)"], color="tab:purple")
        ax3.set_xlabel("Time (UTC)", fontsize=8)
        ax3.set_ylabel("Doppler (kHz)", color="tab:purple", fontsize=8)
        ax3.tick_params(axis="both", labelcolor="tab:purple", labelsize=8)
        ax3.xaxis.set_major_locator(MinuteLocator(interval=1))
        ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax3.grid(True, linestyle=":", alpha=0.7)
        ax3.set_title("Doppler Shift")
        fig.tight_layout()

        canvas_plot_widget = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas_plot_widget.draw()
        canvas_plot_widget.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)

        display_columns = [
            "Time (UTC)",
            "Elevation (°)",
            "Slant Range (km)",
            "Path Loss (dB)",
            "Atmospheric Att (dB)",
            "Pointing Loss (dB)",
            "Off Boresight Angle (°)",
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

    def exit_app(self):
        plt.close("all")
        self.root.quit()
        self.root.destroy()
        sys.exit()

    def show_antenna_pattern(self):
        popup = tk.Toplevel(self.root)
        popup.title("Antenna Pattern")
        popup.geometry("600x400")
        popup.configure(bg=PALETTE["bg"])

        from calculations import ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS, marker="o", linestyle="-")
        ax.set_xlabel("Angle (degrees)", fontsize=8)
        ax.set_ylabel("Gain (dB)", fontsize=8)
        ax.set_title("Antenna Gain Pattern")
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True, linestyle=":", alpha=0.7)
        canvas_popup = FigureCanvasTkAgg(fig, master=popup)
        canvas_popup.draw()
        canvas_popup.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)

    def load_antenna_pattern_file(self):
        """Load antenna gain pattern data from a CSV file."""
        file_path = filedialog.askopenfilename(
            title="Select Antenna Pattern File",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")],
        )
        if not file_path:
            return
        try:
            load_antenna_pattern(file_path)
            messagebox.showinfo("Success", "Antenna pattern loaded successfully.")
            self.set_analysis_stale()
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
            self.current_table_df.to_csv(file_path, index=False)
            messagebox.showinfo("Export Successful", f"Data exported to {file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export CSV: {e}")

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
        self.root.title("Satellite Link Budget Tool")
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

        style.configure(
            "TLabelframe",
            background=PALETTE["surface"],
            bordercolor=PALETTE["border"],
            relief="solid",
            borderwidth=1,
        )
        style.configure(
            "TLabelframe.Label",
            background=PALETTE["surface"],
            foreground=PALETTE["accent_dark"],
            font=bold_font,
        )

        style.configure(
            "TButton",
            background=PALETTE["surface_alt"],
            foreground=PALETTE["text"],
            padding=(12, 7),
            relief="solid",
            borderwidth=1,
            bordercolor=PALETTE["border"],
            font=default_font,
        )
        style.map(
            "TButton",
            background=[("active", PALETTE["accent_light"]), ("disabled", PALETTE["surface_alt"])],
            foreground=[("disabled", PALETTE["text_muted"])],
            bordercolor=[("active", PALETTE["accent"]), ("!active", PALETTE["border"])],
        )

        style.configure(
            "Red.TButton",
            background=PALETTE["danger"],
            foreground="white",
            font=bold_font,
            padding=(12, 7),
            relief="flat",
            borderwidth=0,
        )
        style.map(
            "Red.TButton",
            background=[("active", PALETTE["danger_dark"]), ("disabled", PALETTE["surface_alt"])],
            foreground=[("disabled", PALETTE["text_muted"])],
        )

        style.configure("TNotebook", background=PALETTE["bg"], borderwidth=0)
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
        style.configure(
            "Treeview.Heading",
            background=PALETTE["header_bg"],
            foreground=PALETTE["header_fg"],
            font=bold_font,
            relief="flat",
            padding=(6, 6),
        )
        style.map("Treeview.Heading", background=[("active", PALETTE["accent_dark"])])
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
        badge = tk.Label(
            header,
            text="LB",
            bg=PALETTE["accent"],
            fg="white",
            font=(base_font_family, 13, "bold"),
            width=3,
            padx=4,
            pady=4,
        )
        badge.pack(side=tk.LEFT)
        title_box = ttk.Frame(header, style="Header.TFrame")
        title_box.pack(side=tk.LEFT, padx=(12, 0))
        ttk.Label(title_box, text="Satellite Link Budget Tool", style="Header.TLabel").pack(anchor="w")
        ttk.Label(
            title_box,
            text="Orbit prediction & downlink/uplink budget analysis",
            style="HeaderSub.TLabel",
        ).pack(anchor="w")

        self.root.grid_rowconfigure(0, weight=0)
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        canvas_frame = ttk.Frame(self.root)
        canvas_frame.grid(row=1, column=0, sticky="nsew")
        canvas_frame.grid_rowconfigure(0, weight=1)
        canvas_frame.grid_columnconfigure(0, weight=1)

        canvas = tk.Canvas(canvas_frame, bg=PALETTE["bg"], highlightthickness=0)
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar_y = ttk.Scrollbar(canvas_frame, orient="vertical", command=canvas.yview)
        scrollbar_x = ttk.Scrollbar(canvas_frame, orient="horizontal", command=canvas.xview)
        scrollbar_y.grid(row=0, column=1, sticky="ns")
        scrollbar_x.grid(row=1, column=0, sticky="ew")
        canvas.configure(yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)

        def _on_mousewheel(event):
            if sys.platform == "darwin":
                canvas.yview_scroll(-1 * event.delta, "units")
            else:
                canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _on_shift_mousewheel(event):
            if sys.platform == "darwin":
                canvas.xview_scroll(-1 * event.delta, "units")
            else:
                canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")

        if sys.platform == "darwin":
            canvas.bind_all("<Mousewheel>", _on_mousewheel)
            canvas.bind_all("<Shift-Mousewheel>", _on_shift_mousewheel)
        else:
            canvas.bind_all("<MouseWheel>", _on_mousewheel)
            canvas.bind_all("<Shift-MouseWheel>", _on_shift_mousewheel)

        scrollable_frame = ttk.Frame(canvas)
        canvas_window_id = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw", width=canvas.winfo_width())

        def _on_canvas_configure(event):
            canvas.itemconfig(canvas_window_id, width=event.width)
            canvas.configure(scrollregion=canvas.bbox("all"))

        canvas.bind("<Configure>", _on_canvas_configure)
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        main_frame = ttk.Frame(scrollable_frame, padding=15)
        main_frame.pack(fill="both", expand=True)
        main_frame.grid_columnconfigure(0, weight=1)

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
        ttk.Button(tle_frame, text="Load TLE from file", command=self.load_tle_from_file).grid(
            row=2, column=0, columnspan=2, sticky="w", padx=5, pady=(6, 0)
        )
        tle_frame.grid_columnconfigure(1, weight=1)

        # --- Observation Settings Section ---
        obs_frame = ttk.LabelFrame(main_frame, text="Observation Settings", padding=10)
        obs_frame.pack(fill=tk.X, pady=5)
        ttk.Label(obs_frame, text="Date (YYYY-MM-DD)").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.date_entry = ttk.Entry(obs_frame, width=15)
        self.date_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.date_entry.insert(0, datetime.now(timezone.utc).strftime("%Y-%m-%d"))
        self.date_entry.bind("<KeyRelease>", lambda event: self.set_analysis_stale())
        ttk.Label(obs_frame, text="Ground Station").grid(row=0, column=2, sticky="w", padx=5, pady=2)
        self.gs_var = tk.StringVar(value="")
        self.gs_menu = ttk.Combobox(obs_frame, textvariable=self.gs_var, values=[], state="readonly", width=15)
        self.gs_menu.grid(row=0, column=3, sticky="w", padx=5, pady=2)
        self.gs_menu.bind("<<ComboboxSelected>>", lambda event: self.set_analysis_stale())
        ttk.Button(
            obs_frame,
            text="Load Ground Stations",
            command=self.load_ground_stations_from_file,
        ).grid(row=0, column=4, sticky="w", padx=5, pady=2)
        ttk.Button(obs_frame, text="Load Parameters", command=self.load_parameters_from_file).grid(
            row=0, column=5, sticky="w", padx=5, pady=2
        )
        self.gs_file_var = tk.StringVar(value="Ground stations: none loaded")
        ttk.Label(obs_frame, textvariable=self.gs_file_var, foreground=PALETTE["text_muted"]).grid(
            row=1, column=0, columnspan=6, sticky="w", padx=5, pady=(6, 0)
        )
        self.param_file_var = tk.StringVar(value="Parameters: none loaded")
        ttk.Label(obs_frame, textvariable=self.param_file_var, foreground=PALETTE["text_muted"]).grid(
            row=2, column=0, columnspan=6, sticky="w", padx=5, pady=(2, 0)
        )
        obs_frame.grid_columnconfigure(1, weight=1)
        obs_frame.grid_columnconfigure(3, weight=1)
        obs_frame.grid_columnconfigure(4, weight=1)
        obs_frame.grid_columnconfigure(5, weight=1)

        param_tabs = ttk.Notebook(main_frame)
        param_tabs.pack(fill=tk.BOTH, pady=5, expand=True)

        downlink_tab = ttk.Frame(param_tabs)
        param_tabs.add(downlink_tab, text="Downlink")

        uplink_tab = ttk.Frame(param_tabs)
        param_tabs.add(uplink_tab, text="Uplink")

        param_container = ttk.Frame(downlink_tab)
        param_container.pack(fill=tk.X, expand=False, padx=5, pady=5)
        param_container.grid_columnconfigure(0, weight=1)
        param_container.grid_columnconfigure(1, weight=1)
        param_container.grid_columnconfigure(2, weight=1)
        param_container.grid_columnconfigure(3, weight=1)
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
        atm_frame.grid(row=0, column=2, sticky="nsew", padx=(0, 10), pady=5)
        ttk.Label(atm_frame, text="Link Availability [%]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.LA_entry = ttk.Entry(atm_frame, width=15)
        self.LA_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(atm_frame, text="Other Attenuations [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.other_att_entry = ttk.Entry(atm_frame, width=15)
        self.other_att_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        self.atm_label_var = tk.StringVar(
            value=f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}° El: N/A"
        )
        ttk.Label(atm_frame, textvariable=self.atm_label_var).grid(row=2, column=0, columnspan=2, sticky="w", pady=(10, 0), padx=5)
        atm_frame.grid_columnconfigure(1, weight=1)

        # Baseband Parameters frame
        baseband_frame = ttk.LabelFrame(param_container, text="Baseband Parameters", padding=10)
        baseband_frame.grid(row=0, column=3, sticky="nsew", pady=5)
        ttk.Label(baseband_frame, text="Bit Rate [Mbps]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.bitrate_entry = ttk.Entry(baseband_frame, width=15)
        self.bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Roll-off Factor").grid(row=0, column=2, sticky="w", padx=5, pady=2)
        self.rolloff_entry = ttk.Entry(baseband_frame, width=15)
        self.rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Demodulator Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.demod_loss_entry = ttk.Entry(baseband_frame, width=15)
        self.demod_loss_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(baseband_frame, text="Overhead (Conv. + RS)").grid(row=1, column=2, sticky="w", padx=5, pady=2)
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
        self.bitrate_entry.bind("<KeyRelease>", self.update_link_budget_derived)
        self.rolloff_entry.bind("<KeyRelease>", self.update_link_budget_derived)
        self.overhead_entry.bind("<KeyRelease>", self.update_link_budget_derived)
        self.spectral_eff_entry.bind("<KeyRelease>", self.update_link_budget_derived)
        self.update_link_budget_derived()

        downlink_btn_frame = ttk.LabelFrame(
            downlink_tab, text="Downlink Actions", padding=10
        )
        downlink_btn_frame.pack(fill=tk.X, pady=(0, 10), padx=5)
        self.start_refresh_button = ttk.Button(
            downlink_btn_frame,
            text="Start/Refresh Analysis",
            command=self.run_analysis,
            style="Red.TButton",
        )
        self.start_refresh_button.pack(side=tk.LEFT, padx=5)
        self.recalc_button = ttk.Button(
            downlink_btn_frame,
            text="Calculate DL Link Budget",
            command=self.recalculate_link_budget,
            state="disabled",
        )
        self.recalc_button.pack(side=tk.LEFT, padx=5)
        ttk.Button(downlink_btn_frame, text="Load Antenna Pattern", command=self.load_antenna_pattern_file).pack(
            side=tk.LEFT, padx=5
        )
        ttk.Button(downlink_btn_frame, text="Show Antenna Gain", command=self.show_antenna_pattern).pack(
            side=tk.LEFT, padx=5
        )

        downlink_results_container = ttk.Frame(downlink_tab)
        downlink_results_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))

        btn_frame = ttk.Frame(downlink_results_container)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Export Table CSV", command=self.export_table_csv).pack(
            side=tk.LEFT, padx=5
        )
        ttk.Button(btn_frame, text="Exit", command=self.exit_app).pack(side=tk.RIGHT, padx=5)

        contact_frame = ttk.LabelFrame(
            downlink_results_container, text="Contact Windows (UTC Time)", padding=10
        )
        contact_frame.pack(fill=tk.X, pady=(0, 5))
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

        self.plot_frame = ttk.Frame(downlink_results_container)
        self.plot_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))

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
        ttk.Label(uplink_baseband_frame, text="UL Bit Rate [Mbps]").grid(
            row=0, column=0, sticky="w", padx=5, pady=2
        )
        self.uplink_bitrate_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_baseband_frame, text="UL Roll-off Factor").grid(
            row=0, column=2, sticky="w", padx=5, pady=2
        )
        self.uplink_rolloff_entry = ttk.Entry(uplink_baseband_frame, width=15)
        self.uplink_rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
        ttk.Label(uplink_baseband_frame, text="UL Overhead (Conv. + RS)").grid(
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
            entry.bind("<KeyRelease>", lambda event: self.set_uplink_budget_stale())

        for entry in [
            self.uplink_bitrate_entry,
            self.uplink_rolloff_entry,
            self.uplink_overhead_entry,
            self.uplink_spectral_eff_entry,
        ]:
            entry.bind("<KeyRelease>", self.on_uplink_baseband_change)

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

        uplink_results_frame = ttk.LabelFrame(uplink_tab, text="Uplink Results", padding=10)
        uplink_results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        self.uplink_table_frame = ttk.Frame(uplink_results_frame)
        self.uplink_table_frame.pack(fill=tk.BOTH, expand=True)

        icon_path = os.path.join(base_path, "Satellite.ico")
        if os.path.isfile(icon_path):
            try:
                self.root.iconbitmap(icon_path)
            except tk.TclError:
                pass


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
