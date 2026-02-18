"""Tkinter GUI for the satellite link budget tool."""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
from datetime import datetime, timedelta, timezone
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MinuteLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
from skyfield.api import load, EarthSatellite, wgs84
import astropy.units as u
import os
import sys
import json
import webbrowser

import calculations
from cesium_export import (
    cesiumpy_status,
    export_cesium_bundle,
    preview_cesium_view,
    slice_times,
)

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

# Global variables
contact_windows = []
df_all = pd.DataFrame()
analysis_needs_refresh = True
uplink_needs_refresh = False
current_table_df = pd.DataFrame()
current_gs_file = ""
gs_menu: ttk.Combobox | None = None
gs_file_var: tk.StringVar | None = None
param_file_var: tk.StringVar | None = None
uplink_bitrate_entry: ttk.Entry | None = None
uplink_rolloff_entry: ttk.Entry | None = None
uplink_overhead_entry: ttk.Entry | None = None
uplink_spectral_eff_entry: ttk.Entry | None = None
info_bitrate_ul_var: tk.StringVar | None = None
channel_bw_ul_var: tk.StringVar | None = None
uplink_table_frame: ttk.Frame | None = None
uplink_recalc_button: ttk.Button | None = None


def load_tle_from_file():
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

    tle1_entry.delete(0, tk.END)
    tle1_entry.insert(0, tle1)
    tle2_entry.delete(0, tk.END)
    tle2_entry.insert(0, tle2)
    set_analysis_stale()


def load_parameters_from_file():
    """Load all configurable parameters from a JSON file.

    The JSON structure should provide keys such as ``eirp_sat_dbw`` or
    ``frequency_ghz``. See the README for the full schema and example. Missing
    keys are ignored so the user can provide only the fields they need.
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
        "eirp_sat_dbw": eirp_sat_entry,
        "eirp_gs_dbw": eirp_gs_entry,
        "frequency_ghz": freq_entry,
        "uplink_frequency_ghz": uplink_freq_entry,
        "c_io_dbhz": CIo_entry,
        "gt_gs_dbk": gt_gs_entry,
        "gt_sat_dbk": gt_sat_entry,
        "antenna_diameter_m": d_gs_entry,
        "link_availability_pct": LA_entry,
        "other_attenuations_db": other_att_entry,
        "uplink_other_attenuations_db": other_att_ul_entry,
        "bitrate_mbps": bitrate_entry,
        "rolloff": rolloff_entry,
        "demod_loss_db": demod_loss_entry,
        "demod_loss_sat_db": demod_loss_ul_entry,
        "overhead": overhead_entry,
        "spectral_efficiency_bpshz": spectral_eff_entry,
        "uplink_bitrate_mbps": uplink_bitrate_entry,
        "uplink_rolloff": uplink_rolloff_entry,
        "uplink_overhead": uplink_overhead_entry,
        "uplink_spectral_efficiency_bpshz": uplink_spectral_eff_entry,
    }

    for key, entry in field_map.items():
        if key in payload and entry is not None:
            entry.delete(0, tk.END)
            entry.insert(0, str(payload[key]))

    if param_file_var is not None:
        param_file_var.set(f"Parameters: {os.path.abspath(file_path)}")

    update_link_budget_derived()
    set_analysis_stale()


def load_ground_stations_from_file():
    """Let the user pick a ground station catalogue and refresh the UI."""

    global current_gs_file
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

    current_gs_file = os.path.abspath(file_path)
    if gs_file_var is not None:
        gs_file_var.set(f"Ground stations: {current_gs_file}")
    if gs_menu is not None:
        gs_menu["values"] = list(stations.keys())
    if stations:
        selected = gs_var.get()
        if selected not in stations:
            gs_var.set(next(iter(stations)))
    set_analysis_stale()


def clear_uplink_table():
    """Remove any existing uplink-only table widgets."""

    if uplink_table_frame is None:
        return
    for widget in uplink_table_frame.winfo_children():
        widget.destroy()


def set_analysis_stale():
    """Mark the analysis as stale so the user must recompute."""
    global analysis_needs_refresh
    if not analysis_needs_refresh:
        analysis_needs_refresh = True
        start_refresh_button.config(style="Red.TButton")
        recalc_button.config(state="disabled")
        clear_plot_and_table()
        clear_uplink_table()
        contact_listbox.delete(0, tk.END)


def set_uplink_budget_stale():
    """Mark the uplink-only link budget as stale."""

    global uplink_needs_refresh
    if not uplink_needs_refresh:
        uplink_needs_refresh = True
        clear_uplink_table()
        if uplink_recalc_button is not None:
            uplink_recalc_button.config(style="Red.TButton")


def set_uplink_budget_fresh():
    """Reset the uplink-only link budget stale flag."""

    global uplink_needs_refresh
    if uplink_needs_refresh:
        uplink_needs_refresh = False
        if uplink_recalc_button is not None:
            uplink_recalc_button.config(style="TButton")


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


# --- Main Analysis Function ---

def run_analysis():
    """Run a full orbital and link budget analysis using GUI inputs.

    The function reads TLE data, ground station parameters and link budget
    settings from the GUI fields. It then propagates the orbit for a 24 hour
    period with a 30 second step, computes the link budget for each step using
    :func:`calculate_link_budget_parameters` and stores the results in the
    global :data:`df_all` DataFrame. Detected contact windows are populated in
    :data:`contact_windows` and displayed in the GUI.

    No parameters are accepted; all values are taken from the GUI widgets. The
    function updates the global state and returns ``None``.
    """
    
    global contact_windows, df_all, analysis_needs_refresh, atm_label_var

    tle1 = tle1_entry.get().strip()
    tle2 = tle2_entry.get().strip()
    gs_name = gs_var.get()
    date_str = date_entry.get().strip()

    if not all([tle1, tle2, gs_name, date_str]):
        messagebox.showerror("Error", "All TLE, Ground Station, and Date fields must be filled.")
        return

    try:
        obs_date = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    except ValueError:
        messagebox.showerror("Error", "Invalid date format. Please use YYYY-MM-DD.")
        return

    try:
        freq = float(freq_entry.get()) * u.GHz
        link_availability = float(LA_entry.get())
        p = 100.0 - link_availability
        d_gs = float(d_gs_entry.get())
        other_att = float(other_att_entry.get() or 0.0)
        eirp_sat = float(eirp_sat_entry.get())
        gt_gs = float(gt_gs_entry.get())
        cisat_db = float(CIo_entry.get()) if CIo_entry.get().strip() else None
        cisat_lin = 10 ** (cisat_db / 10.0) if cisat_db is not None else None
        bitrate = float(bitrate_entry.get()) * 1e6
        rolloff = float(rolloff_entry.get())
        demod_loss = float(demod_loss_entry.get())
        overhead = float(overhead_entry.get())
        uplink_freq_val = _get_optional_float(uplink_freq_entry)
        uplink_freq = uplink_freq_val * u.GHz if uplink_freq_val is not None else None
        eirp_gs = _get_optional_float(eirp_gs_entry)
        gt_sat = _get_optional_float(gt_sat_entry)
        demod_loss_ul = _get_optional_float(demod_loss_ul_entry)
        other_att_ul = _get_optional_float(other_att_ul_entry)
        uplink_bitrate_val = _get_optional_float(uplink_bitrate_entry) if uplink_bitrate_entry else None
        uplink_overhead = _get_optional_float(uplink_overhead_entry) if uplink_overhead_entry else None
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
    atm_label_var.set(
        f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}\u00b0 El: {atm_att:.3f}"
    )
    results_list = []
    contact_windows.clear()
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
            contact_windows.append((contact_start_time, contact_end_time, pass_direction))
            in_contact = False
    if in_contact:
        pass_direction = _classify_pass_direction(sat, ts, contact_start_time, times[-1])
        contact_windows.append((contact_start_time, times[-1], pass_direction))

    df_all = pd.DataFrame(results_list)
    display_contacts(contact_windows)
    if contact_windows:
        contact_listbox.selection_set(0)
        on_contact_select(None)
    else:
        clear_plot_and_table()

    analysis_needs_refresh = False
    start_refresh_button.config(style="TButton")
    recalc_button.config(state="!disabled")


# --- Recalculate Link Budget Function ---

def recalculate_link_budget():
    """Recompute the link budget for previously generated contact times.

    This helper function does not re-run the orbital propagation. Instead it
    uses the existing time stamps stored in :data:`df_all` and recalculates the
    link budget with the current parameters from the GUI. The updated results
    replace the contents of :data:`df_all` and the plots/tables shown in the
    interface are refreshed accordingly.

    All configuration values are taken from the GUI widgets and no value is
    returned.
    """
    global df_all, atm_label_var
    if df_all.empty:
        messagebox.showwarning("Warning", "Please generate passes first.")
        return
    if analysis_needs_refresh:
        messagebox.showwarning("Warning", "Please refresh the analysis first.")
        return
    try:
        freq = float(freq_entry.get()) * u.GHz
        link_availability = float(LA_entry.get())
        p = 100.0 - link_availability
        d_gs = float(d_gs_entry.get())
        other_att = float(other_att_entry.get() or 0.0)
        eirp_sat = float(eirp_sat_entry.get())
        gt_gs = float(gt_gs_entry.get())
        cisat_db = float(CIo_entry.get()) if CIo_entry.get().strip() else None
        cisat_lin = 10 ** (cisat_db / 10.0) if cisat_db is not None else None
        bitrate = float(bitrate_entry.get()) * 1e6
        demod_loss = float(demod_loss_entry.get())
        overhead = float(overhead_entry.get())
        uplink_freq_val = _get_optional_float(uplink_freq_entry)
        uplink_freq = uplink_freq_val * u.GHz if uplink_freq_val is not None else None
        eirp_gs = _get_optional_float(eirp_gs_entry)
        gt_sat = _get_optional_float(gt_sat_entry)
        demod_loss_ul = _get_optional_float(demod_loss_ul_entry)
        other_att_ul = _get_optional_float(other_att_ul_entry)
        uplink_bitrate_val = _get_optional_float(uplink_bitrate_entry) if uplink_bitrate_entry else None
        uplink_overhead = _get_optional_float(uplink_overhead_entry) if uplink_overhead_entry else None
    except ValueError as e:
        messagebox.showerror("Input Error", f"Invalid numerical input: {e}.")
        return

    gs_name = gs_var.get()
    lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]
    alt_gs_km = alt_gs_m / 1000.0

    ts = load.timescale()
    tle1 = tle1_entry.get().strip()
    tle2 = tle2_entry.get().strip()
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
    atm_label_var.set(
        f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}\u00b0 El: {atm_att:.3f}"
    )
    updated_results = []
    times_dt = df_all["Time (UTC)"].tolist()
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

    df_all = pd.DataFrame(updated_results)
    selection = contact_listbox.curselection()
    if not selection and len(contact_windows) > 0:
        contact_listbox.selection_set(0)
        selection = (0,)

    if selection:
        on_contact_select(None)
    else:
        clear_plot_and_table()


def clear_plot_and_table():
    for widget in plot_frame.winfo_children():
        widget.destroy()
    for widget in table_frame.winfo_children():
        widget.destroy()


def display_contacts(contacts):
    contact_listbox.delete(0, tk.END)
    filtered_contacts = []
    for contact in contacts:
        if len(contact) == 3:
            start, end, pass_direction = contact
        else:
            start, end = contact
            pass_direction = "A"
        if not df_all.empty:
            mask = (df_all["Time (UTC)"] >= start) & (df_all["Time (UTC)"] <= end)
            df_segment = df_all[mask]
            if (
                not df_segment.empty
                and df_segment["Elevation (°)"].max() >= MIN_ELEVATION_DEG
            ):
                filtered_contacts.append((start, end, pass_direction))
    global contact_windows
    contact_windows = filtered_contacts
    for i, (start, end, pass_direction) in enumerate(contact_windows):
        contact_listbox.insert(
            tk.END,
            f"Contact {i+1}: {start.strftime('%H:%M:%S')} - {end.strftime('%H:%M:%S')} ({pass_direction})",
        )


def on_contact_select(event):
    selection = contact_listbox.curselection()
    if not selection:
        clear_plot_and_table()
        clear_uplink_table()
        return

    idx = selection[0]
    start, end, _ = contact_windows[idx]
    mask = (df_all["Time (UTC)"] >= start) & (df_all["Time (UTC)"] <= end)
    df_pass = df_all[mask].copy()
    global current_table_df
    current_table_df = df_pass

    clear_uplink_table()


    clear_plot_and_table()
    fig, (ax1, ax2, ax3) = plt.subplots(
        1,
        3,
        figsize=(15, 4),
        gridspec_kw={"width_ratios": [1, 1, 1]},
    )
 
    tk_bg_color = root.cget("bg")
    try:
        r, g, b = root.winfo_rgb(tk_bg_color)
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

    canvas_plot_widget = FigureCanvasTkAgg(fig, master=plot_frame)
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
    table = ttk.Treeview(table_frame, columns=display_columns, show="headings")
    for col in display_columns:
        table.heading(col, text=col)
        if col == "Time (UTC)":
            table.column(col, anchor="center", width=120)
        else:
            table.column(col, anchor="center", width=100)
    for _, row in df_pass.iterrows():
        formatted_time = row["Time (UTC)"].strftime("%Y-%m-%d %H:%M:%S")
        values = [formatted_time] + [row.get(col, "") for col in display_columns[1:]]
        table.insert("", "end", values=values)
    table_vscroll = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=table.yview)
    table_hscroll = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=table.xview)
    table.configure(yscrollcommand=table_vscroll.set, xscrollcommand=table_hscroll.set)
    table_vscroll.pack(side=tk.RIGHT, fill=tk.Y)
    table_hscroll.pack(side=tk.BOTTOM, fill=tk.X)
    table.pack(fill=tk.BOTH, expand=True)


def exit_app():
    plt.close("all")
    root.quit()
    root.destroy()
    sys.exit()


def show_antenna_pattern():
    popup = tk.Toplevel(root)
    popup.title("Antenna Pattern")
    popup.geometry("600x400")

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


def load_antenna_pattern_file():
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
        set_analysis_stale()
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load antenna pattern: {e}")


def export_table_csv():
    """Export the currently displayed table to a CSV file."""
    if current_table_df.empty:
        messagebox.showwarning("Warning", "Please select a contact window first.")
        return
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )
    if not file_path:
        return
    try:
        current_table_df.to_csv(file_path, index=False)
        messagebox.showinfo("Export Successful", f"Data exported to {file_path}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to export CSV: {e}")


def export_cesium_view():
    """Export a Cesium CZML/HTML bundle for the selected contact window."""

    if analysis_needs_refresh:
        messagebox.showwarning("Warning", "Please refresh the analysis first.")
        return

    if df_all.empty:
        messagebox.showwarning("Warning", "No analysis data available to export.")
        return

    selection = contact_listbox.curselection()
    if not selection:
        messagebox.showwarning("Warning", "Please select a contact window first.")
        return

    idx = selection[0]
    start, end, _ = contact_windows[idx]
    times = slice_times(df_all["Time (UTC)"].tolist(), start, end)
    if not times:
        messagebox.showwarning("Warning", "No timestamps found for this contact window.")
        return

    output_base = filedialog.asksaveasfilename(
        defaultextension=".html",
        filetypes=[("HTML files", "*.html"), ("All files", "*.*")],
        title="Save Cesium HTML",
    )
    if not output_base:
        return

    tle1 = tle1_entry.get().strip()
    tle2 = tle2_entry.get().strip()
    gs_name = gs_var.get()
    lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]

    try:
        paths = export_cesium_bundle(
            tle1,
            tle2,
            gs_name,
            lat_gs,
            lon_gs,
            alt_gs_m,
            times,
            output_base,
        )
    except Exception as exc:
        messagebox.showerror("Cesium Export", f"Unable to export Cesium bundle: {exc}")
        return

    messagebox.showinfo(
        "Cesium Export",
        "Export completed. Opening HTML in your browser.\n"
        f"HTML: {paths.html_path}\nCZML: {paths.czml_path}",
    )
    webbrowser.open(f"file://{paths.html_path}")


def preview_cesium_view_gui():
    """Render a Cesium preview (requires CesiumPy)."""

    if analysis_needs_refresh:
        messagebox.showwarning("Warning", "Please refresh the analysis first.")
        return

    if df_all.empty:
        messagebox.showwarning("Warning", "No analysis data available to export.")
        return

    selection = contact_listbox.curselection()
    if not selection:
        messagebox.showwarning("Warning", "Please select a contact window first.")
        return

    idx = selection[0]
    start, end, _ = contact_windows[idx]
    times = slice_times(df_all["Time (UTC)"].tolist(), start, end)
    if not times:
        messagebox.showwarning("Warning", "No timestamps found for this contact window.")
        return

    tle1 = tle1_entry.get().strip()
    tle2 = tle2_entry.get().strip()
    gs_name = gs_var.get()
    lat_gs, lon_gs, alt_gs_m = calculations.GROUND_STATIONS[gs_name]

    def _run_preview(google_api_key: str | None = None):
        return preview_cesium_view(
            tle1,
            tle2,
            gs_name,
            lat_gs,
            lon_gs,
            alt_gs_m,
            times,
            google_api_key=google_api_key,
        )

    try:
        paths = _run_preview()
    except Exception as exc:
        if "Google requires" in str(exc) or "google_api_key" in str(exc):
            key = simpledialog.askstring(
                "Cesium Preview",
                "Inserisci la Google Maps API key per CesiumPy:",
                show="*",
            )
            if key:
                try:
                    paths = _run_preview(google_api_key=key)
                except Exception as retry_exc:
                    exc = retry_exc
                else:
                    messagebox.showinfo(
                        "Cesium Preview",
                        "Opening Cesium preview in your browser.\n"
                        f"HTML: {paths.html_path}",
                    )
                    webbrowser.open(f"file://{paths.html_path}")
                    return
        available, detail = cesiumpy_status()
        google_hint = (
            "\nTip: set the CESIUMPY_GOOGLE_API_KEY environment variable "
            "to your Google Maps API key (and verify the key is valid)."
        )
        messagebox.showerror(
            "Cesium Preview",
            "CesiumPy preview is not available.\n"
            f"Details: {exc}\n"
            f"CesiumPy status: {'available' if available else 'missing'} ({detail})\n"
            "Tip: use 'Export Cesium View' to create HTML/CZML without CesiumPy."
            f"{google_hint}",
        )
        return

    messagebox.showinfo(
        "Cesium Preview",
        "Opening Cesium preview in your browser.\n"
        f"HTML: {paths.html_path}",
    )
    webbrowser.open(f"file://{paths.html_path}")


def calculate_ul_link_budget():
    """Render an uplink-only link budget table for the selected contact window."""

    if analysis_needs_refresh:
        messagebox.showwarning("Warning", "Please refresh the analysis first.")
        return

    if uplink_needs_refresh:
        recalculate_link_budget()

    if current_table_df.empty:
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
    missing_columns = required_columns - set(current_table_df.columns)
    if missing_columns:
        messagebox.showerror(
            "Data Error",
            "Uplink link budget values are not available. Please run the analysis again.",
        )
        return

    if uplink_table_frame is None:
        messagebox.showerror("UI Error", "Uplink results area is not initialised.")
        return

    clear_uplink_table()

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

    table = ttk.Treeview(uplink_table_frame, columns=display_columns, show="headings")
    for col in display_columns:
        table.heading(col, text=col)
        if col == "Time (UTC)":
            table.column(col, anchor="center", width=150)
        else:
            table.column(col, anchor="center", width=130)

    for _, row in current_table_df.iterrows():
        formatted_time = row["Time (UTC)"].strftime("%Y-%m-%d %H:%M:%S")
        values = [formatted_time] + [_format_cell(row.get(col, "")) for col in display_columns[1:]]
        table.insert("", "end", values=values)

    table_vscroll = ttk.Scrollbar(uplink_table_frame, orient=tk.VERTICAL, command=table.yview)
    table_hscroll = ttk.Scrollbar(uplink_table_frame, orient=tk.HORIZONTAL, command=table.xview)
    table.configure(yscrollcommand=table_vscroll.set, xscrollcommand=table_hscroll.set)
    table_vscroll.pack(side=tk.RIGHT, fill=tk.Y)
    table_hscroll.pack(side=tk.BOTTOM, fill=tk.X)
    table.pack(fill=tk.BOTH, expand=True)
    set_uplink_budget_fresh()


def _update_baseband_info(
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


def update_link_budget_derived(*args):
    _update_baseband_info(
        bitrate_entry,
        rolloff_entry,
        overhead_entry,
        spectral_eff_entry,
        info_bitrate_var,
        channel_bw_var,
        "",
    )

    if info_bitrate_ul_var is not None and channel_bw_ul_var is not None and all(
        entry is not None
        for entry in (
            uplink_bitrate_entry,
            uplink_rolloff_entry,
            uplink_overhead_entry,
            uplink_spectral_eff_entry,
        )
    ):
        _update_baseband_info(
            uplink_bitrate_entry,
            uplink_rolloff_entry,
            uplink_overhead_entry,
            uplink_spectral_eff_entry,
            info_bitrate_ul_var,
            channel_bw_ul_var,
            "UL",
        )


def on_uplink_baseband_change(event=None):
    """Refresh derived UL metrics and mark uplink budget for recalculation."""

    update_link_budget_derived()
    set_uplink_budget_stale()


# --- GUI Setup ---

def setup_gui():
    """Initialise and launch the Tkinter graphical interface."""
    global root, tle1_entry, tle2_entry, date_entry, gs_var, freq_entry
    global LA_entry, d_gs_entry, other_att_entry, eirp_sat_entry
    global eirp_gs_entry, uplink_freq_entry, gt_gs_entry, gt_sat_entry
    global demod_loss_ul_entry, other_att_ul_entry
    global CIo_entry, bitrate_entry, rolloff_entry, demod_loss_entry
    global overhead_entry, spectral_eff_entry, atm_label_var, info_bitrate_var, channel_bw_var
    global uplink_bitrate_entry, uplink_rolloff_entry, uplink_overhead_entry, uplink_spectral_eff_entry
    global info_bitrate_ul_var, channel_bw_ul_var
    global recalc_button, contact_listbox, plot_frame, table_frame
    global start_refresh_button, gs_menu, gs_file_var, param_file_var
    global uplink_table_frame, uplink_recalc_button

    root = tk.Tk()
    root.title("Satellite Link Budget Tool")
    root.geometry("1200x950")

    style = ttk.Style()
    style.theme_use("clam")
    style.configure("TButton", padding=(10, 6), relief="raised", borderwidth=2)
    style.map(
        "TButton",
        relief=[("pressed", "sunken"), ("!pressed", "raised")],
    )
    style.configure(
        "Red.TButton",
        background="red",
        foreground="white",
        font=("TkDefaultFont", 10, "bold"),
        padding=(10, 6),
        relief="raised",
        borderwidth=2,
    )
    style.map(
        "Red.TButton",
        background=[("active", "darkred"), ("!active", "red")],
        foreground=[("active", "white"), ("!active", "white")],
        relief=[("pressed", "sunken"), ("!pressed", "raised")],
    )

    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)

    canvas_frame = ttk.Frame(root)
    canvas_frame.grid(row=0, column=0, sticky="nsew")
    canvas_frame.grid_rowconfigure(0, weight=1)
    canvas_frame.grid_columnconfigure(0, weight=1)

    canvas = tk.Canvas(canvas_frame)
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
    tle1_entry = ttk.Entry(tle_frame, width=80)
    tle1_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    tle1_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    ttk.Label(tle_frame, text="TLE Line 2").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    tle2_entry = ttk.Entry(tle_frame, width=80)
    tle2_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    tle2_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    ttk.Button(tle_frame, text="Load TLE from file", command=load_tle_from_file).grid(
        row=2, column=0, columnspan=2, sticky="w", padx=5, pady=(6, 0)
    )
    tle_frame.grid_columnconfigure(1, weight=1)

    # --- Observation Settings Section ---
    obs_frame = ttk.LabelFrame(main_frame, text="Observation Settings", padding=10)
    obs_frame.pack(fill=tk.X, pady=5)
    ttk.Label(obs_frame, text="Date (YYYY-MM-DD)").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    date_entry = ttk.Entry(obs_frame, width=15)
    date_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
    date_entry.insert(0, datetime.now(timezone.utc).strftime("%Y-%m-%d"))
    date_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    ttk.Label(obs_frame, text="Ground Station").grid(row=0, column=2, sticky="w", padx=5, pady=2)
    gs_var = tk.StringVar(value="")
    gs_menu = ttk.Combobox(obs_frame, textvariable=gs_var, values=[], state="readonly", width=15)
    gs_menu.grid(row=0, column=3, sticky="w", padx=5, pady=2)
    gs_menu.bind("<<ComboboxSelected>>", lambda event: set_analysis_stale())
    ttk.Button(
        obs_frame,
        text="Load Ground Stations",
        command=load_ground_stations_from_file,
    ).grid(row=0, column=4, sticky="w", padx=5, pady=2)
    ttk.Button(obs_frame, text="Load Parameters", command=load_parameters_from_file).grid(
        row=0, column=5, sticky="w", padx=5, pady=2
    )
    gs_file_var = tk.StringVar(value="Ground stations: none loaded")
    ttk.Label(obs_frame, textvariable=gs_file_var, foreground="gray25").grid(
        row=1, column=0, columnspan=6, sticky="w", padx=5, pady=(6, 0)
    )
    param_file_var = tk.StringVar(value="Parameters: none loaded")
    ttk.Label(obs_frame, textvariable=param_file_var, foreground="gray25").grid(
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
    eirp_sat_entry = ttk.Entry(sat_frame, width=15)
    eirp_sat_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(sat_frame, text="Frequency [GHz]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    freq_entry = ttk.Entry(sat_frame, width=15)
    freq_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(sat_frame, text="C/Io [dBHz]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
    CIo_entry = ttk.Entry(sat_frame, width=15)
    CIo_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
    sat_frame.grid_columnconfigure(1, weight=1)

    # Ground Station Parameters frame
    gs_frame = ttk.LabelFrame(param_container, text="Ground Station Parameters", padding=10)
    gs_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 10), pady=5)
    ttk.Label(gs_frame, text="G/T GS [dB/K]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    gt_gs_entry = ttk.Entry(gs_frame, width=15)
    gt_gs_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(gs_frame, text="Antenna Diameter GS [m]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    d_gs_entry = ttk.Entry(gs_frame, width=15)
    d_gs_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    gs_frame.grid_columnconfigure(1, weight=1)

    # Atmospheric Attenuations frame
    atm_frame = ttk.LabelFrame(param_container, text="Atmospheric Attenuations", padding=10)
    atm_frame.grid(row=0, column=2, sticky="nsew", padx=(0, 10), pady=5)
    ttk.Label(atm_frame, text="Link Availability [%]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    LA_entry = ttk.Entry(atm_frame, width=15)
    LA_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(atm_frame, text="Other Attenuations [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    other_att_entry = ttk.Entry(atm_frame, width=15)
    other_att_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    atm_label_var = tk.StringVar(
        value=f"Atmospheric Att (dB) @ {MIN_ELEVATION_DEG:g}° El: N/A"
    )
    ttk.Label(atm_frame, textvariable=atm_label_var).grid(row=2, column=0, columnspan=2, sticky="w", pady=(10, 0), padx=5)
    atm_frame.grid_columnconfigure(1, weight=1)

    # Baseband Parameters frame
    baseband_frame = ttk.LabelFrame(param_container, text="Baseband Parameters", padding=10)
    baseband_frame.grid(row=0, column=3, sticky="nsew", pady=5)
    ttk.Label(baseband_frame, text="Bit Rate [Mbps]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    bitrate_entry = ttk.Entry(baseband_frame, width=15)
    bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Roll-off Factor").grid(row=0, column=2, sticky="w", padx=5, pady=2)
    rolloff_entry = ttk.Entry(baseband_frame, width=15)
    rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Demodulator Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    demod_loss_entry = ttk.Entry(baseband_frame, width=15)
    demod_loss_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Overhead (Conv. + RS)").grid(row=1, column=2, sticky="w", padx=5, pady=2)
    overhead_entry = ttk.Entry(baseband_frame, width=15)
    overhead_entry.grid(row=1, column=3, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Spectral Efficiency [b/s/Hz]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
    spectral_eff_entry = ttk.Entry(baseband_frame, width=15)
    spectral_eff_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
    info_bitrate_var = tk.StringVar(value="Info Bit Rate [Mbps]: N/A")
    channel_bw_var = tk.StringVar(value="Channel BW [MHz]: N/A")
    ttk.Label(baseband_frame, textvariable=info_bitrate_var).grid(row=3, column=0, columnspan=2, sticky="w", pady=(5, 0), padx=5)
    ttk.Label(baseband_frame, textvariable=channel_bw_var).grid(row=3, column=2, columnspan=2, sticky="w", padx=5)
    baseband_frame.grid_columnconfigure(1, weight=1)
    baseband_frame.grid_columnconfigure(3, weight=1)
    bitrate_entry.bind("<KeyRelease>", update_link_budget_derived)
    rolloff_entry.bind("<KeyRelease>", update_link_budget_derived)
    overhead_entry.bind("<KeyRelease>", update_link_budget_derived)
    spectral_eff_entry.bind("<KeyRelease>", update_link_budget_derived)
    update_link_budget_derived()

    downlink_btn_frame = ttk.LabelFrame(
        downlink_tab, text="Downlink Actions", padding=10
    )
    downlink_btn_frame.pack(fill=tk.X, pady=(0, 10), padx=5)
    start_refresh_button = ttk.Button(
        downlink_btn_frame,
        text="Start/Refresh Analysis",
        command=run_analysis,
        style="Red.TButton",
    )
    start_refresh_button.pack(side=tk.LEFT, padx=5)
    recalc_button = ttk.Button(
        downlink_btn_frame,
        text="Calculate DL Link Budget",
        command=recalculate_link_budget,
        state="disabled",
    )
    recalc_button.pack(side=tk.LEFT, padx=5)
    ttk.Button(downlink_btn_frame, text="Load Antenna Pattern", command=load_antenna_pattern_file).pack(
        side=tk.LEFT, padx=5
    )
    ttk.Button(downlink_btn_frame, text="Show Antenna Gain", command=show_antenna_pattern).pack(
        side=tk.LEFT, padx=5
    )

    downlink_results_container = ttk.Frame(downlink_tab)
    downlink_results_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))

    btn_frame = ttk.Frame(downlink_results_container)
    btn_frame.pack(fill=tk.X, pady=(0, 10))
    ttk.Button(btn_frame, text="Export Table CSV", command=export_table_csv).pack(
        side=tk.LEFT, padx=5
    )
    ttk.Button(btn_frame, text="Export Cesium View", command=export_cesium_view).pack(
        side=tk.LEFT, padx=5
    )
    ttk.Button(
        btn_frame,
        text="Preview Cesium (CesiumPy)",
        command=preview_cesium_view_gui,
    ).pack(side=tk.LEFT, padx=5)
    ttk.Button(btn_frame, text="Exit", command=exit_app).pack(side=tk.RIGHT, padx=5)

    contact_frame = ttk.LabelFrame(
        downlink_results_container, text="Contact Windows (UTC Time)", padding=10
    )
    contact_frame.pack(fill=tk.X, pady=(0, 5))
    contact_listbox = tk.Listbox(contact_frame, height=6, exportselection=False)
    contact_listbox.pack(fill=tk.X, expand=True)
    contact_listbox.bind("<<ListboxSelect>>", on_contact_select)

    plot_frame = ttk.Frame(downlink_results_container)
    plot_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))

    downlink_results_frame = ttk.LabelFrame(
        downlink_results_container, text="Downlink Results", padding=10
    )
    downlink_results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))
    table_frame = ttk.Frame(downlink_results_frame)
    table_frame.pack(fill=tk.BOTH, expand=True)

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
    eirp_gs_entry = ttk.Entry(uplink_gs_frame, width=15)
    eirp_gs_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(uplink_gs_frame, text="Uplink Frequency [GHz]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    uplink_freq_entry = ttk.Entry(uplink_gs_frame, width=15)
    uplink_freq_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    uplink_gs_frame.grid_columnconfigure(1, weight=1)

    uplink_sat_frame = ttk.LabelFrame(uplink_container, text="Satellite Parameters", padding=10)
    uplink_sat_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 10), pady=5)
    ttk.Label(uplink_sat_frame, text="G/T SAT [dB/K]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    gt_sat_entry = ttk.Entry(uplink_sat_frame, width=15)
    gt_sat_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(uplink_sat_frame, text="Sat Demod Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    demod_loss_ul_entry = ttk.Entry(uplink_sat_frame, width=15)
    demod_loss_ul_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    uplink_sat_frame.grid_columnconfigure(1, weight=1)

    uplink_atm_frame = ttk.LabelFrame(uplink_container, text="Uplink Attenuations", padding=10)
    uplink_atm_frame.grid(row=0, column=2, sticky="nsew", pady=5)
    ttk.Label(uplink_atm_frame, text="Uplink Other Attenuations [dB]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    other_att_ul_entry = ttk.Entry(uplink_atm_frame, width=15)
    other_att_ul_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    uplink_atm_frame.grid_columnconfigure(1, weight=1)

    uplink_baseband_frame = ttk.LabelFrame(
        uplink_container, text="Uplink Baseband Parameters", padding=10
    )
    uplink_baseband_frame.grid(row=1, column=0, columnspan=3, sticky="nsew", pady=(0, 5))
    ttk.Label(uplink_baseband_frame, text="UL Bit Rate [Mbps]").grid(
        row=0, column=0, sticky="w", padx=5, pady=2
    )
    uplink_bitrate_entry = ttk.Entry(uplink_baseband_frame, width=15)
    uplink_bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(uplink_baseband_frame, text="UL Roll-off Factor").grid(
        row=0, column=2, sticky="w", padx=5, pady=2
    )
    uplink_rolloff_entry = ttk.Entry(uplink_baseband_frame, width=15)
    uplink_rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
    ttk.Label(uplink_baseband_frame, text="UL Overhead (Conv. + RS)").grid(
        row=1, column=0, sticky="w", padx=5, pady=2
    )
    uplink_overhead_entry = ttk.Entry(uplink_baseband_frame, width=15)
    uplink_overhead_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(uplink_baseband_frame, text="UL Spectral Efficiency [b/s/Hz]").grid(
        row=1, column=2, sticky="w", padx=5, pady=2
    )
    uplink_spectral_eff_entry = ttk.Entry(uplink_baseband_frame, width=15)
    uplink_spectral_eff_entry.grid(row=1, column=3, sticky="ew", padx=5, pady=2)
    info_bitrate_ul_var = tk.StringVar(value="UL Info Bit Rate [Mbps]: N/A")
    channel_bw_ul_var = tk.StringVar(value="UL Channel BW [MHz]: N/A")
    ttk.Label(uplink_baseband_frame, textvariable=info_bitrate_ul_var).grid(
        row=2, column=0, columnspan=2, sticky="w", pady=(5, 0), padx=5
    )
    ttk.Label(uplink_baseband_frame, textvariable=channel_bw_ul_var).grid(
        row=2, column=2, columnspan=2, sticky="w", padx=5
    )

    uplink_baseband_frame.grid_columnconfigure(1, weight=1)
    uplink_baseband_frame.grid_columnconfigure(3, weight=1)

    for entry in [
        eirp_gs_entry,
        uplink_freq_entry,
        gt_sat_entry,
        demod_loss_ul_entry,
        other_att_ul_entry,
    ]:
        entry.bind("<KeyRelease>", lambda event: set_uplink_budget_stale())

    for entry in [
        uplink_bitrate_entry,
        uplink_rolloff_entry,
        uplink_overhead_entry,
        uplink_spectral_eff_entry,
    ]:
        entry.bind("<KeyRelease>", on_uplink_baseband_change)

    uplink_actions_frame = ttk.LabelFrame(
        uplink_tab, text="Uplink Actions", padding=10
    )
    uplink_actions_frame.pack(fill=tk.X, pady=10, padx=5)
    uplink_recalc_button = ttk.Button(
        uplink_actions_frame,
        text="Calculate UL Link Budget",
        command=calculate_ul_link_budget,
    )
    uplink_recalc_button.pack(side=tk.LEFT, padx=5)

    uplink_results_frame = ttk.LabelFrame(uplink_tab, text="Uplink Results", padding=10)
    uplink_results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
    uplink_table_frame = ttk.Frame(uplink_results_frame)
    uplink_table_frame.pack(fill=tk.BOTH, expand=True)

    root.mainloop()


if __name__ == "__main__":
    setup_gui()
    root.iconbitmap(os.path.join(base_path, "Satellite.ico"))
