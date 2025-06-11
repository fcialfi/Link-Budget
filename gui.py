"""Tkinter GUI for the satellite link budget tool."""
import tkinter as tk
from tkinter import ttk, messagebox
from datetime import datetime, timedelta, timezone
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MinuteLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import sys

from skyfield.api import load, EarthSatellite, wgs84
import astropy.units as u

from .calculations import (
    GROUND_STATIONS,
    calculate_link_budget_parameters,
)

# Global variables
contact_windows = []
df_all = pd.DataFrame()
analysis_needs_refresh = True


def set_analysis_stale():
    """Mark the analysis as stale so the user must recompute."""
    global analysis_needs_refresh
    if not analysis_needs_refresh:
        analysis_needs_refresh = True
        start_refresh_button.config(style="Red.TButton")
        recalc_button.config(state="disabled")
        clear_plot_and_table()
        contact_listbox.delete(0, tk.END)


# --- Main Analysis Function ---

def run_analysis():
    global contact_windows, df_all, analysis_needs_refresh

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
        r001 = float(r001_entry.get())
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
    except ValueError as e:
        messagebox.showerror("Input Error", f"Invalid numerical input: {e}.")
        return

    lat_gs, lon_gs, alt_gs_m = GROUND_STATIONS[gs_name]
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

    diff_at_all_times = (sat - gs).at(sky_times)
    altitudes, azimuths, distances = diff_at_all_times.altaz()

    results_list = []
    contact_windows.clear()
    in_contact = False
    contact_start_time = None

    for i, t_utc_dt in enumerate(times):
        elev = altitudes.degrees[i]
        t_sky = sky_times[i]
        params = calculate_link_budget_parameters(
            t_sky,
            sat,
            gs,
            freq,
            p,
            r001,
            d_gs,
            alt_gs_km,
            eirp_sat,
            gt_gs,
            demod_loss,
            bitrate,
            overhead,
            cisat_lin,
            other_att,
        )
        for key in [
            "Elevation (°)",
            "Slant Range (km)",
            "Path Loss (dB)",
            "Pointing Loss (dB)",
            "Off Boresight Angle (°)",
            "Rx Power (dBW)",
            "C/(No+Io) (dBHz)",
            "Eb/No (dB)",
        ]:
            if params[key] is not None:
                params[key] = round(params[key], 2)
        results_list.append(params)

        is_visible_now = elev >= 5
        if is_visible_now and not in_contact:
            contact_start_time = t_utc_dt
            in_contact = True
        elif not is_visible_now and in_contact:
            contact_end_time = t_utc_dt
            contact_windows.append((contact_start_time, contact_end_time))
            in_contact = False
    if in_contact:
        contact_windows.append((contact_start_time, times[-1]))

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
    global df_all
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
        r001 = float(r001_entry.get())
        d_gs = float(d_gs_entry.get())
        other_att = float(other_att_entry.get() or 0.0)
        eirp_sat = float(eirp_sat_entry.get())
        gt_gs = float(gt_gs_entry.get())
        cisat_db = float(CIo_entry.get()) if CIo_entry.get().strip() else None
        cisat_lin = 10 ** (cisat_db / 10.0) if cisat_db is not None else None
        bitrate = float(bitrate_entry.get()) * 1e6
        demod_loss = float(demod_loss_entry.get())
        overhead = float(overhead_entry.get())
    except ValueError as e:
        messagebox.showerror("Input Error", f"Invalid numerical input: {e}.")
        return

    gs_name = gs_var.get()
    lat_gs, lon_gs, alt_gs_m = GROUND_STATIONS[gs_name]
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

    updated_results = []
    for _, row in df_all.iterrows():
        t_utc_dt = row["Time (UTC)"]
        t_sky = ts.utc(
            t_utc_dt.year,
            t_utc_dt.month,
            t_utc_dt.day,
            t_utc_dt.hour,
            t_utc_dt.minute,
            t_utc_dt.second,
        )
        params = calculate_link_budget_parameters(
            t_sky,
            sat,
            gs,
            freq,
            p,
            r001,
            d_gs,
            alt_gs_km,
            eirp_sat,
            gt_gs,
            demod_loss,
            bitrate,
            overhead,
            cisat_lin,
            other_att,
        )
        for key in [
            "Elevation (°)",
            "Slant Range (km)",
            "Path Loss (dB)",
            "Pointing Loss (dB)",
            "Off Boresight Angle (°)",
            "Rx Power (dBW)",
            "C/(No+Io) (dBHz)",
            "Eb/No (dB)",
        ]:
            if params[key] is not None:
                params[key] = round(params[key], 2)
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
    for start, end in contacts:
        if not df_all.empty:
            mask = (df_all["Time (UTC)"] >= start) & (df_all["Time (UTC)"] <= end)
            df_segment = df_all[mask]
            if not df_segment.empty and df_segment["Elevation (°)"].max() >= 5:
                filtered_contacts.append((start, end))
    global contact_windows
    contact_windows = filtered_contacts
    for i, (start, end) in enumerate(contact_windows):
        contact_listbox.insert(
            tk.END,
            f"Contact {i+1}: {start.strftime('%H:%M:%S')} - {end.strftime('%H:%M:%S')}",
        )


def on_contact_select(event):
    selection = contact_listbox.curselection()
    if not selection:
        clear_plot_and_table()
        return

    idx = selection[0]
    start, end = contact_windows[idx]
    mask = (df_all["Time (UTC)"] >= start) & (df_all["Time (UTC)"] <= end)
    df_pass = df_all[mask].copy()

    clear_plot_and_table()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1]})
    tk_bg_color = root.cget("bg")
    try:
        r, g, b = root.winfo_rgb(tk_bg_color)
        mpl_bg_color = f"#{r // 256:02x}{g // 256:02x}{b // 256:02x}"
    except tk.TclError:
        mpl_bg_color = "#F0F0F0"
    fig.patch.set_facecolor(mpl_bg_color)

    ax1.plot(df_pass["Time (UTC)"], df_pass["Elevation (°)"], label="Elevation (°)", color="tab:blue")
    ax1.set_ylabel("Elevation (°)", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.set_xlabel("Time (UTC)")
    ax1.xaxis.set_major_locator(MinuteLocator(interval=1))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.grid(True, linestyle=":", alpha=0.7)
    ax1.set_title("Elevation Angle")
    fig.autofmt_xdate(rotation=30)

    ax2_1 = ax2
    ax2_2 = ax2.twinx()
    ax2_1.plot(df_pass["Time (UTC)"], df_pass["Eb/No (dB)"], label="Eb/No (dB)", color="tab:red")
    ax2_1.set_ylabel("Eb/No (dB)", color="tab:red")
    ax2_1.tick_params(axis="y", labelcolor="tab:red")
    ax2_2.plot(df_pass["Time (UTC)"], df_pass["C/(No+Io) (dBHz)"], label="C/(N+Io) (dBHz)", color="tab:green", linestyle="--")
    ax2_2.set_ylabel("C/(No+Io) (dBHz)", color="tab:green")
    ax2_2.tick_params(axis="y", labelcolor="tab:green")
    ax2.set_xlabel("Time (UTC)")
    ax2.xaxis.set_major_locator(MinuteLocator(interval=1))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax2.grid(True, linestyle=":", alpha=0.7)
    ax2.set_title("Link Quality (Eb/No & C/(No+Io))")
    lines, labels = ax2_1.get_legend_handles_labels()
    lines2, labels2 = ax2_2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="best")
    fig.tight_layout()

    canvas_plot_widget = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas_plot_widget.draw()
    canvas_plot_widget.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)

    display_columns = [
        "Time (UTC)",
        "Elevation (°)",
        "Slant Range (km)",
        "Path Loss (dB)",
        "Pointing Loss (dB)",
        "Off Boresight Angle (°)",
        "Rx Power (dBW)",
        "C/(No+Io) (dBHz)",
        "Eb/No (dB)",
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

    from .calculations import ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS, marker="o", linestyle="-")
    ax.set_xlabel("Angle (degrees)")
    ax.set_ylabel("Gain (dB)")
    ax.set_title("Antenna Gain Pattern")
    ax.grid(True, linestyle=":", alpha=0.7)
    canvas_popup = FigureCanvasTkAgg(fig, master=popup)
    canvas_popup.draw()
    canvas_popup.get_tk_widget().pack(side=tk.LEFT, anchor="nw", fill=tk.BOTH, expand=True)


def update_link_budget_derived(*args):
    try:
        bit_rate_mbps = float(bitrate_entry.get())
        roll_off = float(rolloff_entry.get())
        overhead = float(overhead_entry.get())
        info_bit_rate_mbps = bit_rate_mbps / overhead if overhead != 0 else 0
        channel_bw_mhz = bit_rate_mbps * (1 + roll_off) / 2
        info_bitrate_var.set(f"Info Bit Rate [Mbps]: {info_bit_rate_mbps:.3f}")
        channel_bw_var.set(f"Channel BW [MHz] (QPSK): {channel_bw_mhz:.3f}")
    except ValueError:
        info_bitrate_var.set("Info Bit Rate [Mbps]: N/A")
        channel_bw_var.set("Channel BW [MHz] (QPSK): N/A")


# --- GUI Setup ---

def setup_gui():
    global root, tle1_entry, tle2_entry, date_entry, gs_var, freq_entry
    global LA_entry, r001_entry, d_gs_entry, other_att_entry, eirp_sat_entry
    global gt_gs_entry, CIo_entry, bitrate_entry, rolloff_entry, demod_loss_entry
    global overhead_entry, atm_label_var, info_bitrate_var, channel_bw_var
    global recalc_button, contact_listbox, plot_frame, table_frame
    global start_refresh_button

    root = tk.Tk()
    root.title("Satellite Link Budget Tool")
    root.geometry("1200x950")

    style = ttk.Style()
    style.theme_use("clam")
    style.configure("Red.TButton", background="red", foreground="white", font=("TkDefaultFont", 10, "bold"))
    style.map(
        "Red.TButton",
        background=[("active", "darkred"), ("!active", "red")],
        foreground=[("active", "white"), ("!active", "white")],
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
    tle1_entry.insert(0, "1 60543U 24149CD  25156.09113434  .00001223  00000-0 13313-3 0  9994")
    tle1_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    ttk.Label(tle_frame, text="TLE Line 2").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    tle2_entry = ttk.Entry(tle_frame, width=80)
    tle2_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    tle2_entry.insert(0, "2 60543  97.7181 231.4010 0000622 289.0774  71.0376 14.89798875 43525")
    tle2_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    tle_frame.grid_columnconfigure(1, weight=1)

    # --- Observation Settings Section ---
    obs_frame = ttk.LabelFrame(main_frame, text="Observation Settings", padding=10)
    obs_frame.pack(fill=tk.X, pady=5)
    ttk.Label(obs_frame, text="Date (YYYY-MM-DD)").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    date_entry = ttk.Entry(obs_frame, width=15)
    date_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
    date_entry.insert(0, "2025-06-16")
    date_entry.bind("<KeyRelease>", lambda event: set_analysis_stale())
    ttk.Label(obs_frame, text="Ground Station").grid(row=0, column=2, sticky="w", padx=5, pady=2)
    gs_var = tk.StringVar(value="Darmstadt")
    gs_menu = ttk.Combobox(obs_frame, textvariable=gs_var, values=list(GROUND_STATIONS.keys()), state="readonly", width=15)
    gs_menu.grid(row=0, column=3, sticky="w", padx=5, pady=2)
    gs_menu.bind("<<ComboboxSelected>>", lambda event: set_analysis_stale())
    obs_frame.grid_columnconfigure(1, weight=1)
    obs_frame.grid_columnconfigure(3, weight=1)

    param_container = ttk.Frame(main_frame)
    param_container.pack(fill=tk.X, pady=5, expand=True)
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
    eirp_sat_entry.insert(0, "11")
    eirp_sat_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(sat_frame, text="Frequency [GHz]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    freq_entry = ttk.Entry(sat_frame, width=15)
    freq_entry.insert(0, "1.707")
    freq_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(sat_frame, text="C/Io [dBHz]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
    CIo_entry = ttk.Entry(sat_frame, width=15)
    CIo_entry.insert(0, "")
    CIo_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
    sat_frame.grid_columnconfigure(1, weight=1)

    # Ground Station Parameters frame
    gs_frame = ttk.LabelFrame(param_container, text="Ground Station Parameters", padding=10)
    gs_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 10), pady=5)
    ttk.Label(gs_frame, text="G/T GS [dB/K]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    gt_gs_entry = ttk.Entry(gs_frame, width=15)
    gt_gs_entry.insert(0, "8.5")
    gt_gs_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(gs_frame, text="Antenna Diameter GS [m]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    d_gs_entry = ttk.Entry(gs_frame, width=15)
    d_gs_entry.insert(0, "3.7")
    d_gs_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    gs_frame.grid_columnconfigure(1, weight=1)

    # Atmospheric Attenuations frame
    atm_frame = ttk.LabelFrame(param_container, text="Atmospheric Attenuations", padding=10)
    atm_frame.grid(row=0, column=2, sticky="nsew", padx=(0, 10), pady=5)
    ttk.Label(atm_frame, text="Link Availability [%]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    LA_entry = ttk.Entry(atm_frame, width=15)
    LA_entry.insert(0, "99")
    LA_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(atm_frame, text="R001 [mm/h]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    r001_entry = ttk.Entry(atm_frame, width=15)
    r001_entry.insert(0, "12.5")
    r001_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(atm_frame, text="Other Attenuations [dB]").grid(row=2, column=0, sticky="w", padx=5, pady=2)
    other_att_entry = ttk.Entry(atm_frame, width=15)
    other_att_entry.insert(0, "0.0")
    other_att_entry.grid(row=2, column=1, sticky="ew", padx=5, pady=2)
    atm_label_var = tk.StringVar(value="Atmospheric Att (dB) @ 5° El: N/A")
    ttk.Label(atm_frame, textvariable=atm_label_var).grid(row=3, column=0, columnspan=2, sticky="w", pady=(10, 0), padx=5)
    atm_frame.grid_columnconfigure(1, weight=1)

    # Baseband Parameters frame
    baseband_frame = ttk.LabelFrame(param_container, text="Baseband Parameters", padding=10)
    baseband_frame.grid(row=0, column=3, sticky="nsew", pady=5)
    ttk.Label(baseband_frame, text="Bit Rate [Mbps]").grid(row=0, column=0, sticky="w", padx=5, pady=2)
    bitrate_entry = ttk.Entry(baseband_frame, width=15)
    bitrate_entry.insert(0, "3.57")
    bitrate_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Roll-off Factor").grid(row=0, column=2, sticky="w", padx=5, pady=2)
    rolloff_entry = ttk.Entry(baseband_frame, width=15)
    rolloff_entry.insert(0, "0.45")
    rolloff_entry.grid(row=0, column=3, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Demodulator Loss [dB]").grid(row=1, column=0, sticky="w", padx=5, pady=2)
    demod_loss_entry = ttk.Entry(baseband_frame, width=15)
    demod_loss_entry.insert(0, "1.5")
    demod_loss_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=2)
    ttk.Label(baseband_frame, text="Overhead (Conv. + RS)").grid(row=1, column=2, sticky="w", padx=5, pady=2)
    overhead_entry = ttk.Entry(baseband_frame, width=15)
    overhead_entry.insert(0, "2.29")
    overhead_entry.grid(row=1, column=3, sticky="ew", padx=5, pady=2)
    info_bitrate_var = tk.StringVar(value="Info Bit Rate [Mbps]: N/A")
    channel_bw_var = tk.StringVar(value="Channel BW [MHz] (QPSK): N/A")
    ttk.Label(baseband_frame, textvariable=info_bitrate_var).grid(row=2, column=0, columnspan=2, sticky="w", pady=(5, 0), padx=5)
    ttk.Label(baseband_frame, textvariable=channel_bw_var).grid(row=2, column=2, columnspan=2, sticky="w", padx=5)
    baseband_frame.grid_columnconfigure(1, weight=1)
    baseband_frame.grid_columnconfigure(3, weight=1)
    bitrate_entry.bind("<KeyRelease>", update_link_budget_derived)
    rolloff_entry.bind("<KeyRelease>", update_link_budget_derived)
    overhead_entry.bind("<KeyRelease>", update_link_budget_derived)
    update_link_budget_derived()

    # --- Buttons frame ---
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(fill=tk.X, pady=10)
    start_refresh_button = ttk.Button(btn_frame, text="Start/Refresh Analysis", command=run_analysis, style="Red.TButton")
    start_refresh_button.pack(side=tk.LEFT, padx=5)
    recalc_button = ttk.Button(btn_frame, text="Recalculate Link Budget", command=recalculate_link_budget, state="disabled")
    recalc_button.pack(side=tk.LEFT, padx=5)
    ttk.Button(btn_frame, text="Show Antenna Gain", command=show_antenna_pattern).pack(side=tk.LEFT, padx=5)
    ttk.Button(btn_frame, text="Exit", command=exit_app).pack(side=tk.RIGHT, padx=5)

    contact_frame = ttk.LabelFrame(main_frame, text="Contact Windows (UTC Time)", padding=10)
    contact_frame.pack(fill=tk.X, pady=5)
    contact_listbox = tk.Listbox(contact_frame, height=6, exportselection=False)
    contact_listbox.pack(fill=tk.X, expand=True)
    contact_listbox.bind("<<ListboxSelect>>", on_contact_select)

    plot_frame = ttk.Frame(main_frame)
    plot_frame.pack(fill=tk.BOTH, expand=True, pady=5)

    table_frame = ttk.Frame(main_frame)
    table_frame.pack(fill=tk.BOTH, expand=True, pady=5)

    root.mainloop()


if __name__ == "__main__":
    setup_gui()
