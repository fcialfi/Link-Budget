"""Chapters of the user manual that document every parameter and how it is computed.

Called by build_manual.py, which passes its own module (``m``) for the
building blocks (H1, H2, P, TABLE, ...) so the chapter and figure numbers
stay in one sequence. The worked examples are computed here, at build time,
with the functions of the program, so the numbers printed in the manual are
the ones the program gives.
"""

import datetime
import re
import warnings

import numpy as np


def _fmt(value, digits=2):
    return f"{value:.{digits}f}"


# --------------------------------------------------------------------------
# Worked examples (computed with the program's own functions)
# --------------------------------------------------------------------------
# An invented mission and ground station, used by every example of the manual.
EXAMPLE = {
    "station": "Example station",
    "lat": 64.0, "lon": 20.0, "alt_m": 300.0, "d_gs": 3.7,
    "freq_ghz": 2.25, "elevation_deg": 5.0, "altitude_km": 650.0, "availability_pct": 99.9,
    "tx_power_w": 4.0, "circuit_loss_db": 0.8, "vswr": 1.3,
    # Satellite antenna pattern: angle from nadir [deg], gain [dBi].
    "pattern": ([0, 20, 40, 50, 60, 65, 70], [4.0, 3.6, 2.4, 1.4, 0.2, -0.5, -1.5]),
    "coded_mbps": 4.0, "interleaving": 4, "rolloff": 0.35, "spectral_eff": 2.0,
    "required_ebno_db": 2.5, "required_margin_db": 3.0,
    "demod_loss_db": 1.5, "modulation_degradation_db": 0.3,
    "s4_ref": 0.2, "tx_axial_ratio_db": 2.0, "rx_axial_ratio_db": 1.0, "multipath_db": 0.2,
    "rx_depointing_deg": 0.1,
    "gt_dbk": 13.5, "gt_ref_elevation_deg": 10.0, "tsys_k": 160.0, "radome_db": 0.3,
    "pfd_low": -154.0, "pfd_high": -144.0,
    "tolerances": [
        {"name": "Transmitter Power", "key": "tx_power", "favourable_db": 0.5, "adverse_db": 0.5,
         "distribution": "TRI"},
        {"name": "Atmospheric / Propagation Loss", "key": "propagation", "favourable_db": 1.0,
         "adverse_db": 1.0, "distribution": "GAU"},
    ],
}


def _static_example():
    """The invented downlink, computed as the Static LB · Downlink tab computes it."""
    import astropy.units as u

    import calculations as c

    e = EXAMPLE
    c.set_antenna_pattern(np.array(e["pattern"][0], dtype=float), np.array(e["pattern"][1], dtype=float))
    eta = c.off_nadir_angle_deg(e["elevation_deg"], e["altitude_km"], e["alt_m"] / 1000.0)
    gain = c.satellite_antenna_gain_dbi(eta)
    framing = c.ccsds_tm_overheads(e["interleaving"], True, False)
    rate = e["coded_mbps"] * 1e6
    bandwidth = rate * (1 + e["rolloff"]) / e["spectral_eff"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        budget = c.calculate_fixed_elevation_link_budget(
            freq=e["freq_ghz"] * u.GHz, elevation_deg=e["elevation_deg"], sat_altitude_km=e["altitude_km"],
            lat_gs=e["lat"], lon_gs=e["lon"], alt_gs_km=e["alt_m"] / 1000.0, d_gs=e["d_gs"], eirp=None,
            gt=e["gt_dbk"], demod_loss=e["demod_loss_db"], bitrate=rate, overhead=framing["total"],
            other_att=0.0, pointing_loss_db=0.0, link_availability_pct=e["availability_pct"],
            include_scintillation=False, required_ebno=e["required_ebno_db"], tx_power_w=e["tx_power_w"],
            antenna_circuit_loss_db=e["circuit_loss_db"], vswr=e["vswr"], antenna_gain_dbi=gain,
            ionospheric_s4_ref=e["s4_ref"], polarisation_loss_db=None, multipath_loss_db=e["multipath_db"],
            modulation_degradation_db=e["modulation_degradation_db"], occupied_bandwidth_hz=bandwidth,
            pfd_limit_dbw_m2_4khz=e["pfd_low"], pfd_limit_high_dbw_m2_4khz=e["pfd_high"],
            tx_axial_ratio_db=e["tx_axial_ratio_db"], rx_axial_ratio_db=e["rx_axial_ratio_db"],
            rx_depointing_deg=e["rx_depointing_deg"], rx_beamwidth_3db_deg=None,
            formatting_overhead=framing["formatting"], tolerances=e["tolerances"],
            gt_model={"ref_elevation_deg": e["gt_ref_elevation_deg"], "tsys_ref_k": e["tsys_k"],
                      "radome_loss_db": e["radome_db"]},
        )
    return eta, gain, bandwidth, budget, framing


def _gt_example():
    """G/T of the invented station versus elevation, clear sky and rain, as the program computes it."""
    import calculations as c

    e = EXAMPLE
    p = 100.0 - e["availability_pct"]
    alt_km = e["alt_m"] / 1000.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sky_ref = c.gt_reference_sky_noise_k(e["lat"], e["lon"], e["freq_ghz"], e["d_gs"], alt_km,
                                             e["gt_ref_elevation_deg"])
        sky_ref_zenith = c.gt_reference_sky_noise_k(e["lat"], e["lon"], e["freq_ghz"], e["d_gs"], alt_km, 90.0)
        rows = []
        for el in (2.0, 5.0, 10.0, 20.0, 45.0, 90.0):
            clear = c.atmospheric_attenuation_contributions(
                e["lat"], e["lon"], e["freq_ghz"], el, p, e["d_gs"], alt_km,
                include_rain=False, include_clouds=False, include_scintillation=False)
            rain = c.atmospheric_attenuation_contributions(
                e["lat"], e["lon"], e["freq_ghz"], el, p, e["d_gs"], alt_km, include_scintillation=False)
            a_clear = clear["gas"]
            a_rain = rain["gas"] + rain["cloud"] + rain["rain"]
            t_clear, t_rain = c.sky_noise_temperature_k(a_clear), c.sky_noise_temperature_k(a_rain)
            gt_clear = c.ground_station_gt_dbk(e["gt_dbk"], e["tsys_k"], sky_ref, t_clear, e["radome_db"])
            gt_rain = c.ground_station_gt_dbk(e["gt_dbk"], e["tsys_k"], sky_ref, t_rain, e["radome_db"])
            gt_zenith_ref = c.ground_station_gt_dbk(e["gt_dbk"], e["tsys_k"], sky_ref_zenith, t_clear,
                                                    e["radome_db"])
            rows.append((el, a_clear, t_clear, e["tsys_k"] - sky_ref + t_clear, gt_clear, a_rain, t_rain,
                         gt_rain, gt_zenith_ref))
    return sky_ref, sky_ref_zenith, rows


def _statistics_example():
    """One satellite, two 10-minute contacts in 4 hours: the memory and latency model step by step."""
    import pass_statistics as ps

    step = 30.0
    n = 480  # 4 hours
    visible = np.zeros(n, dtype=bool)
    visible[120:140] = True  # 60 to 70 min
    visible[320:340] = True  # 160 to 170 min
    start = datetime.datetime(2026, 1, 1, tzinfo=datetime.timezone.utc)
    day = ps.DaySamples(start, np.where(visible, 30.0, -10.0), visible, np.where(visible, 5.0, np.nan))
    runs = {}
    for at_end in (False, True):
        runs[at_end] = ps.aggregate([day], step, info_bitrate_bps=1.5e6, min_margin_db=3.0,
                                    generated_bps=60e3, latency_requirement=(75.0, 80.0),
                                    latency_at_pass_end=at_end)
    return runs


# --------------------------------------------------------------------------
# Chapters
# --------------------------------------------------------------------------
def chapters(m, root=None):
    """The reference chapters, as a list of flowables."""
    # Section references are written as tokens (§gt, §atm, ...) and numbered here, so they follow the
    # position of these chapters in the manual.
    ref_ch = m._chapter[0] + 1
    sections = ["const", "symbols", "orbit", "stations", "geom", "atm", "gt", "sun", "dl", "ul", "columns", "static",
                "static_results", "derived", "example"]
    tokens = {name: f"{ref_ch}.{k}" for k, name in enumerate(sections, start=1)}
    tokens.update(ref=str(ref_ch), lat=f"{ref_ch + 1}.3", calc=str(m.CALC_CHAPTER[0]), ref_stats=str(ref_ch + 1))

    def T(text):
        for name in sorted(tokens, key=len, reverse=True):
            text = text.replace("§" + name, tokens[name])
        return text

    H1, H2, H3 = m.H1, m.H2, m.H3

    def P(text):
        return m.P(T(text))

    def NOTE(text):
        return m.NOTE(T(text))

    def BUL(items):
        return m.BUL([T(item) for item in items])

    def TABLE(rows, *args, **kwargs):
        return m.TABLE([[T(cell) if isinstance(cell, str) else cell for cell in row] for row in rows],
                       *args, **kwargs)

    def M(tex):
        """A formula in a table cell."""
        return m.MathFlowable(tex, size=9.5)

    def EQ(*rows):
        return m.EQN([(tex, T(note)) for tex, note in rows])

    W4 = [3.0, 4.6, 1.7, 7.2]  # parameter tables: name, key / unit, empty value, use
    W3 = [3.8, 6.4, 6.3]  # result tables: name, formula, notes

    story = []

    # ======================================================================
    story += H1("Parameter and calculation reference")
    story.append(P(
        "This chapter documents every input of the program and every value it computes: where the input is "
        "entered, its key in <b>parameters.json</b>, its unit, what the program assumes when it is left empty, "
        "and the exact formula through which it reaches the results. It is written from the source code; the "
        "Python function that implements each step is named in italics (<i>module.function</i>), so that every "
        "statement can be checked against the program. Chapter §calc is the short summary of the same methods."))
    story.append(P(
        "Notation: log is the base-10 logarithm; quantities in dB, dBW, dBi, dB/K, dBHz add and subtract; "
        "x<sub>lin</sub> = 10<super>x/10</super> converts a dB value to a linear ratio. Angles are in degrees "
        "unless written in radians. 'Empty' means a field left blank or a null value in parameters.json."))

    story += H2("Constants and units")
    story += TABLE([
        ["Quantity", "Value", "Used for"],
        ["Earth radius R", "6378.137 km (WGS-84 equatorial)", "Static geometry, synthetic TLE (altitude = a − R)."],
        ["Earth gravitational parameter μ", "398600.4418 km³/s²", "Synthetic TLE, Doppler rate."],
        ["J2", "1.08262668·10⁻³", "Sun-synchronous inclination, nodal period."],
        ["Sun-synchronous node rate", "2π / (365.2421897 · 86400) rad/s", "Sun-synchronous inclination."],
        ["Earth rotation rate", "7.2921159·10⁻⁵ rad/s", "Doppler rate (station acceleration)."],
        ["Speed of light c", "299 792.458 km/s", "Doppler, wavelength λ = c / f."],
        ["Boltzmann constant k", "1.380649·10⁻²³ J/K; −10 log k = 228.6 dBW/K/Hz", "C/No = Rx power + G/T + 228.6."],
        ["Mean radiating temperature Tm", "275 K (ITU-R P.618 §3)", "Sky noise temperature."],
        ["Cosmic background Tc", "2.7 K", "Sky noise temperature."],
        ["Sun angular radius", "0.2665° → solid angle Ω<sub>sun</sub> = π (0.2665° in rad)² = 6.80·10⁻⁵ sr",
         "Sun noise."],
        ["Geomagnetic north pole", "80.65° N, 72.68° W (IGRF-13, epoch 2020)", "Geomagnetic latitude, "
                                                                              "ionospheric zone."],
        ["GPS L1", "1575.42 MHz", "Default reference frequency of S4."],
        ["Time step", "30 s", "Pass analysis and statistics; days start at 00:00 UTC."],
        ["Data units", "kbps = 10³ bit/s, Mbps = 10⁶ bit/s; MB = 10⁶ bytes = 8·10⁶ bit; GB = 10⁹ bytes",
         "On-board memory in GB (memory_gb × 8·10⁹ bit), volumes in MB."],
        ["Atmospheric models", "ITU-R P.618-13 and the recommendations it calls (P.676 gas, P.840 cloud, "
                               "P.837/P.838/P.839 rain), library itur 0.4", "Atmospheric attenuation."],
    ], [4.0, 6.0, 6.5], bold_cols=(0,))

    # ----------------------------------------------------------------------
    story += H2("List of symbols")
    story.append(P(
        "Every symbol used in the formulas of chapters §calc, §ref and §ref_stats, with its meaning, its unit "
        "and the section where it is defined or computed. Subscripts: GS ground station, SAT satellite, UL "
        "uplink, ref reference, req required."))
    groups = [
        ("Orbit and geometry", [
            (r"R", "Earth equatorial radius, 6378.137", "km", "§const"),
            (r"\mu,\ J_2", "Earth gravitational parameter, second zonal harmonic", "km³/s², –", "§const"),
            (r"a,\ e,\ i,\ \omega", "Semi-major axis, eccentricity, inclination, argument of perigee", "km, –, °, °", "§orbit"),
            (r"h", "Satellite altitude (a − R)", "km", "§orbit"),
            (r"n,\ p", "Mean motion, semi-latus rectum a(1 − e²)", "rad/s, km", "§orbit"),
            (r"\Omega,\ \dot{\Omega}_{\mathrm{SSO}}", "Right ascension of the ascending node (RAAN); node rate of a sun-synchronous orbit", "°, rad/s", "§orbit"),
            (r"M_0,\ \dot{M},\ \dot{\omega}", "Mean anomaly at the epoch; J2 rates of mean anomaly and perigee", "°, rad/s", "§orbit"),
            (r"T", "Nodal period", "min", "§orbit"),
            (r"u,\ \Delta u", "Argument of latitude at the epoch; phasing between planes", "°", "§orbit"),
            (r"N,\ p,\ j", "Satellites per plane; plane index; satellite index in the plane", "–", "§orbit"),
            (r"t_0,\ \mathrm{LTAN}", "Epoch (00:00 UTC of the day); local time of the ascending node", "–, h", "§orbit"),
            (r"\mathrm{el},\ \mathrm{az}", "Elevation and azimuth of the satellite seen from the station", "°", "§stations"),
            (r"\mathrm{el}_{\min},\ \mathrm{el}_{\mathrm{req}},\ m", "Min Elevation; required elevation; horizon mask", "°", "§stations"),
            (r"h_{\mathrm{GS}}", "Station altitude", "km", "§geom"),
            (r"d", "Slant range", "km (m in the PFD)", "§geom"),
            (r"\eta", "Angle from the satellite nadir to the station (off-boresight angle)", "°", "§geom"),
            (r"f_D,\ \dot{f}_D", "Doppler shift and Doppler rate", "Hz, Hz/s", "§geom"),
            (r"\dot{r},\ \ddot{r}", "Range rate and range acceleration", "km/s, km/s²", "§geom"),
            (r"\mathbf{v}_{\mathrm{rel}},\ \hat{\mathbf{r}}", "Relative velocity; unit vector from the station to the satellite", "km/s, –", "§geom"),
            (r"\mathbf{a}_{\mathrm{sat}},\ \mathbf{a}_{\mathrm{GS}}", "Gravity acceleration of the satellite; centripetal acceleration of the station", "km/s²", "§geom"),
            (r"\theta", "Angle between two directions (Sun separation, depointing)", "°", "§geom, §sun"),
        ]),
        ("Propagation and losses", [
            (r"f,\ \lambda", "Carrier frequency; wavelength c / f", "GHz, m", "§dl"),
            (r"p", "Percentage of time the attenuation is exceeded: 100 − availability", "%", "§atm"),
            (r"A,\ A_g,\ A_c,\ A_r,\ A_s", "Total, gas, cloud, rain, tropospheric scintillation attenuation", "dB", "§atm"),
            (r"L_{\mathrm{FS}}", "Free space loss", "dB", "§columns"),
            (r"L_{\mathrm{iono}},\ S_4,\ P_{\mathrm{pp}}", "Ionospheric loss; S4 scintillation index; peak-to-peak fluctuation", "dB, –, dB", "§derived"),
            (r"f_{\mathrm{ref}},\ S_{4,\mathrm{ref}}", "Frequency at which S4 was measured; S4 at that frequency", "MHz, –", "§derived"),
            (r"\phi,\ \lambda,\ \phi_m,\ \phi_p,\ \lambda_p", "Station latitude and longitude; geomagnetic latitude; geomagnetic pole", "°", "§derived"),
            (r"L_{\mathrm{pol}},\ AR,\ r", "Polarisation mismatch loss; axial ratio; voltage axial ratio 10<super>AR/20</super>", "dB, dB, –", "§derived"),
            (r"L_{\mathrm{mp}},\ L_{\mathrm{other}}", "Multipath losses; Other Attenuations", "dB", "§static"),
            (r"L_{\mathrm{link}}", "Static losses applied along the passes: L_iono + L_pol + L_mp", "dB", "§columns"),
            (r"L_{\mathrm{prop}}", "Total propagation loss of the static budget", "dB", "§static_results"),
            (r"L_{\mathrm{point}},\ \theta_p,\ \theta_{3\mathrm{dB}}", "Pointing loss; depointing angle; 3 dB beamwidth 70 λ/D", "dB, °, °", "§static_results"),
        ]),
        ("Antennas, G/T and Sun", [
            (r"P_{\mathrm{Tx}}", "Transmitter power (10 log of the power in W)", "dBW", "§static_results"),
            (r"L_{\mathrm{circuit}},\ L_{\mathrm{VSWR}},\ \Gamma", "Antenna circuit loss; mismatch loss; reflection coefficient", "dB, dB, –", "§static_results"),
            (r"T_{\mathrm{chain}}", "Tx chain without antenna gain: P_Tx − L_circuit − L_VSWR", "dBW", "§columns"),
            (r"G_{\mathrm{ant}},\ G_{\mathrm{eff}}", "Tx antenna gain; effective gain G_ant − L_circuit − L_VSWR", "dBi", "§static_results"),
            (r"G(\eta),\ G_{\mathrm{peak}}", "Satellite antenna pattern gain at η; highest gain of the pattern", "dBi", "§columns"),
            (r"\mathrm{EIRP}", "Effective isotropic radiated power", "dBW", "§columns"),
            (r"D,\ D_{\mathrm{ref}}", "Ground antenna diameter; diameter the G/T GS refers to", "m", "§gt"),
            (r"(G/T)_{\mathrm{GS}},\ (G/T)_{\mathrm{ref}}", "G/T GS entered; G/T of the station at the reference elevation", "dB/K", "§gt"),
            (r"\mathrm{el}_{\mathrm{ref}}", "G/T reference elevation", "°", "§gt"),
            (r"T_{\mathrm{sys}},\ T_{\mathrm{sys,ref}}", "System noise temperature; at the reference elevation, clear sky", "K", "§gt"),
            (r"T_{\mathrm{sky}},\ T_{\mathrm{sky,ref}}", "Sky noise temperature; at the reference elevation, clear sky", "K", "§gt"),
            (r"T_m,\ T_c", "Mean radiating temperature 275 K; cosmic background 2.7 K", "K", "§const"),
            (r"L_{\mathrm{radome}}", "Radome loss", "dB", "§gt"),
            (r"T_b,\ S,\ \Omega_{\odot}", "Sun brightness temperature; solar flux; solar disc solid angle", "K, sfu, sr", "§sun"),
            (r"G_{\mathrm{main}},\ G_{\mathrm{side}}", "Ground antenna main-lobe and side-lobe gain towards the Sun", "dBi", "§sun"),
            (r"\Delta T,\ \Delta_{\mathrm{Sun}}", "Noise added by the Sun; Sun noise (drop of G/T)", "K, dB", "§sun"),
            (r"k", "Boltzmann constant (−228.6 dBW/K/Hz); in §orbit the J2 coefficient", "–", "§const"),
            (r"n,\ L,\ g,\ \lambda_{\odot},\ \varepsilon", "Days from J2000; Sun mean longitude, mean anomaly, ecliptic longitude; obliquity", "d, °", "§sun"),
            (r"\alpha,\ \delta,\ H,\ \mathrm{GMST}", "Sun right ascension, declination, hour angle; sidereal time (α is the roll-off elsewhere)", "°", "§sun"),
            (r"\lambda_{\mathrm{GS}}", "Station longitude", "°", "§sun"),
        ]),
        ("Link budget and data", [
            (r"P_{\mathrm{Rx}}", "Received power", "dBW", "§columns"),
            (r"S/N_0,\ C/N_0,\ C/I_0", "Signal and carrier to noise density; carrier to interference density", "dBHz", "§columns"),
            (r"L_{\mathrm{demod}},\ L_{\mathrm{mod}},\ D_{\mathrm{mod}}", "Demodulator loss; modulation degradation (L_mod = D_mod)", "dB", "§dl, §static"),
            (r"E_b/N_0,\ (E_b/N_0)_{\mathrm{req}}", "Energy per bit to noise density; required", "dB", "§columns"),
            (r"R_c,\ R_{\mathrm{frame}},\ R_{\mathrm{info}}", "Coded, transfer frame and information bit rate", "bit/s", "§dl, §derived"),
            (r"OH,\ F,\ c", "Overhead R_c/R_info; formatting factor R_frame/R_info; coding factor R_c/R_frame", "–", "§derived"),
            (r"I,\ L_{\mathrm{frame}},\ L_{\mathrm{cb}},\ L_{\mathrm{df}},\ L_{\mathrm{user}}", "RS interleaving depth; frame, codeblock, data field and user data lengths", "–, bytes", "§derived"),
            (r"L_{\mathrm{SH}},\ L_{\mathrm{pkt}},\ L_{\mathrm{CUC}},\ L_{\mathrm{PEC}}", "Frame secondary header; packet length; packet data field header; packet error control", "bytes", "§derived"),
            (r"B,\ \alpha,\ \eta_s", "Occupied bandwidth; roll-off; spectral efficiency", "Hz, –, b/s/Hz", "§dl"),
            (r"\mathrm{PFD},\ \mathrm{PFD}_{4\mathrm{kHz}}", "Power flux density; in the 4 kHz reference bandwidth", "dBW/m², dBW/m²/4kHz", "§static_results"),
            (r"L_{\mathrm{low}},\ L_{\mathrm{high}},\ M_{\mathrm{PFD}}", "PFD limits 0–5° and 25–90°; PFD margin", "dBW/m²/4kHz, dB", "§derived"),
        ]),
        ("Margins and tolerances", [
            (r"M,\ M_{\mathrm{nom}}", "Data Recovery Margin E_b/N_0 − (E_b/N_0)_req (nominal)", "dB", "§static_results"),
            (r"a,\ b,\ \mathrm{adv},\ \mathrm{fav}", "Adverse (negative) and favourable limits of a tolerance", "dB", "§derived"),
            (r"\mu_k,\ \sigma_k^2,\ \sigma", "Mean and variance of tolerance k; standard deviation of the margin", "dB, dB², dB", "§derived"),
            (r"M_{\mathrm{mean}},\ M_{-3\sigma},\ M_{\mathrm{RSS}}", "Margin mean, MEAN − 3σ, Worst Case RSS", "dB", "§derived"),
            (r"M_{\mathrm{fav}},\ M_{\mathrm{adv}}", "Margin with all tolerances favourable / adverse", "dB", "§derived"),
        ]),
        ("Statistics", [
            (r"s,\ g,\ t,\ \Delta t", "Satellite, station, sample; time step 30 s", "–, –, –, s", "§ref_stats.2"),
            (r"m(s,g,t),\ M_{\mathrm{req}}", "Margin over the required Eb/No; Required Margin", "dB", "§ref_stats.2"),
            (r"C,\ R_d,\ R_{\mathrm{RT}}", "Capacity of a sample; downlink rate of the stored data; real-time TM rate", "bit, bit/s, bit/s", "§ref_stats.2"),
            (r"G", "Data generated on board: instruments + HKTM + other", "bit/s", "§ref_stats.3"),
            (r"Q,\ Q_{\max}", "Data on board; on-board memory", "bit", "§ref_stats.3"),
            (r"x,\ D", "Bits downloaded in a sample; bits downloaded before it", "bit", "§ref_stats.3"),
            (r"t_{\mathrm{arr}}", "Time the chunk reaches the ground (end of sample or of pass)", "s", "§ref_stats.3"),
            (r"\ell_{\min},\ \ell_{\max},\ \bar{\ell}", "Latency range of a chunk; mean latency", "s, min", "§ref_stats.3"),
            (r"s(L),\ L_X", "Share of the data within L minutes; latency at X %", "–, min", "§ref_stats.3"),
        ]),
    ]
    for title, items in groups:
        story.append(H3(title))
        story += TABLE([["Symbol", "Meaning", "Unit", "Section"]] +
                       [[M(tex), re.sub(r"([A-Za-z0-9)])_([A-Za-z0-9]+)", r"\1<sub>\2</sub>", meaning), unit, where]
                        for tex, meaning, unit, where in items],
                       [3.6, 8.4, 2.6, 1.9], bold_cols=())

    story += H2("Orbit, TLE and constellation")
    story.append(P(
        "<b>Launched satellite</b>: the two TLE lines of <b>tle.txt</b> in the mission folder, or downloaded by "
        "NORAD ID (key <font face='Courier'>norad_id</font>; sources SatNOGS, then N2YO). The orbit is propagated "
        "with SGP4 (Skyfield) at every 30 s sample. The station is a WGS-84 geodetic point at its latitude, "
        "longitude and altitude."))
    story.append(P(
        "<b>Mission not launched yet</b> (custom orbit, <font face='Courier'>\"custom_orbit\": true</font> in "
        "satellites.json): the program writes a synthetic TLE from the <font face='Courier'>orbit</font> section "
        "of parameters.json (<i>tle_source.synthetic_tle</i>):"))
    story += TABLE([
        ["Parameter", "Key / unit", "Empty", "Use"],
        ["Mean altitude", "orbit.altitude_km [km]", "628", "Semi-major axis a = R + altitude (100–40000 km)."],
        ["Eccentricity", "orbit.eccentricity", "0", "e, 0 ≤ e &lt; 0.1."],
        ["Argument of perigee", "orbit.arg_perigee_deg [deg]", "0", "ω (90 for a frozen orbit)."],
        ["Inclination", "orbit.inclination_deg [deg]", "SSO", "i; empty gives the sun-synchronous inclination "
                                                              "below."],
        ["Planes (LTAN)", "orbit.ltan [hh:mm, list]", "22:20", "Local time of the ascending node of each plane; "
                                                             "sets the RAAN."],
        ["Satellites per plane", "orbit.satellites_per_plane", "1", "N satellites evenly spaced in each plane."],
        ["Phasing between planes", "orbit.phasing_deg [deg]", "0", "Shift of the argument of latitude from one "
                                                                    "plane to the next."],
        ["Extra satellites", "orbit.extra_satellites [LTAN@deg, list]", "none", "Single satellites added, e.g. "
                                                                               "19:30@120."],
    ], W4, mono_cols=(1,))
    story.append(P("The elements of the synthetic TLE are computed as follows:"))
    story.append(EQ(
        (r"n = \sqrt{\mu / a^3}, \quad p = a\,(1 - e^2)", "Mean motion and semi-latus rectum."),
        (r"\cos i_{\mathrm{SSO}} = -\dfrac{\dot{\Omega}_{\mathrm{SSO}}}{1.5\, n\, J_2\, (R/p)^2}",
         "Sun-synchronous inclination, when the inclination is empty."),
        (r"\Omega = \alpha_{\odot}(t_0) + 15^{\circ}\,(\mathrm{LTAN} - 12\,\mathrm{h})",
         "RAAN; α<sub>⊙</sub>: right ascension of the Sun at the epoch t<sub>0</sub>."),
        (r"M_0 = u - \omega", "Mean anomaly; u: argument of latitude of the satellite at the epoch."),
        (r"k = \dfrac{3}{4}\, J_2 \left(\dfrac{R}{p}\right)^2", ""),
        (r"\dot{M} = n \left[1 + k \sqrt{1 - e^2}\,(3\cos^2 i - 1)\right], \quad \dot{\omega} = n\,k\,(5\cos^2 i - 1)",
         "J2 secular rates of the mean anomaly and of the perigee."),
        (r"T = \dfrac{2\pi}{\dot{M} + \dot{\omega}}", "Nodal period."),
    ))
    story.append(P(
        "The epoch is 00:00 UTC of the day analysed (of the first day for the statistics); the drag term is 0. "
        "The mean motion written in the TLE is adjusted (up to 5 iterations, to 10⁻¹⁰) until SGP4's own "
        "Ṁ + ω̇ equals 2π / T: SGP4 then flies the J2 nodal period, and the pass times stay right over the "
        "analysed days (a Keplerian mean motion would drift by about a minute per orbit). The Sun right "
        "ascension uses the low-precision formulae of the Astronomical Almanac (section §sun)."))
    story.append(P(
        "<b>Constellation</b> (<i>link_core.constellation_members</i>): with P planes in the LTAN list, the "
        "satellite j = 1…N of plane p = 0…P−1 has the argument of latitude at the epoch"))
    story.append(EQ(
        (r"u_{p,j} = (j - 1)\,\dfrac{360^{\circ}}{N} + p\,\Delta u \quad (\mathrm{mod}\ 360^{\circ})",
         "Δu: phasing between planes."),
    ))
    story.append(P(
        "An extra satellite \"LTAN@u\" is added with that plane and that u, numbered after the regular "
        "satellites of its plane. The satellites are named \"&lt;mission&gt; &lt;LTAN&gt; #j\" and get the "
        "catalogue numbers 99900, 99901, …"))
    try:
        import tle_source as ts

        alt, ecc = 650.0, 0.001
        inc = ts.sso_inclination_deg(alt, ecc)
        period = ts.nodal_period_min(ts.EARTH_RADIUS_KM + alt, ecc, inc)
        story.append(NOTE(
            f"<b>Example</b> (an invented mission): 650 km, e = 0.001, inclination empty. The sun-synchronous "
            f"inclination is {inc:.3f}° and the nodal period {period:.2f} min. With planes 10:30 and 22:30, "
            "N = 3 satellites per plane and phasing 60°, u is 0°, 120°, 240° in the 10:30 plane and 60°, 180°, "
            "300° in the 22:30 plane; an extra satellite \"10:30@90\" becomes \"&lt;mission&gt; 10:30 #4\" "
            "at u = 90°."))
    except Exception:  # pragma: no cover - the example is only illustrative
        pass

    # ----------------------------------------------------------------------
    story += H2("Ground stations and visibility")
    story.append(P(
        "The stations come from <b>ground_stations.csv</b> (<i>calculations._read_ground_stations</i>). Columns:"))
    story += TABLE([
        ["Column", "Unit", "Empty", "Use"],
        ["name, latitude_deg, longitude_deg, altitude_m", "deg, deg, m", "required",
         "Station position: geometry (WGS-84), ITU-R climate maps, Sun position, geomagnetic latitude."],
        ["antenna_diameter_m", "m", "field Antenna Diameter GS, else 3 m",
         "D: ITU-R scintillation, G/T scaling, 3 dB beamwidth, Sun noise antenna model."],
        ["gt_dbk", "dB/K", "mission G/T",
         "The station's own measured G/T (clear sky); replaces the mission G/T GS and its diameter scaling."],
        ["gt_elevation_deg", "deg", "mission G/T Reference Elevation (5° when that is empty too)",
         "Elevation at which gt_dbk was measured."],
        ["tsys_k", "K", "mission System Noise Temperature", "System noise temperature belonging to gt_dbk."],
        ["antennas", "", "1", "Antennas of the site, backups included. The statistics use 1 per site unless "
                              "Use backup antennas is ticked."],
    ], [4.6, 1.6, 3.4, 6.9], mono_cols=(0,))
    story.append(P(
        "<b>Minimum elevation</b> (Observation Settings, key <font face='Courier'>min_elevation_deg</font>, empty: "
        "5°, allowed −5° to 90°) and the optional <b>horizon mask</b> "
        "(<font face='Courier'>horizon_masks/&lt;station&gt;.txt</font> or .csv: azimuth from North clockwise "
        "and minimum elevation, one point per line, interpolated linearly with wrap-around at North) give the "
        "visibility of each sample:"))
    story.append(EQ(
        (r"\mathrm{el}_{\mathrm{req}}(t) = \max\left[\mathrm{el}_{\min},\ m(\mathrm{az}(t))\right]",
         "el<sub>min</sub>: Min Elevation; m: horizon mask (0 without one)."),
        (r"\mathrm{contact}(t) \Leftrightarrow \mathrm{el}(t) \geq \mathrm{el}_{\mathrm{req}}(t)", ""),
    ))
    story.append(P(
        "A <b>contact</b> is a run of consecutive samples in contact. In the pass list it starts at its first "
        "sample in contact and ends at the first sample out of contact (or at the last sample of the day); "
        "<b>A</b>/<b>D</b> (ascending/descending) compares the sub-satellite latitude at the start and the end. In "
        "the statistics a contact of k samples lasts k × 30 s."))

    # ----------------------------------------------------------------------
    story += H2("Geometry along the passes")
    story.append(P(
        "At every sample Skyfield gives the topocentric position of the satellite seen from the station "
        "(<i>calculations.prepare_topocentric_data</i>):"))
    story += TABLE([
        ["Output", "Formula", "Notes"],
        ["Elevation, azimuth, slant range d", "Topocentric altitude, azimuth and distance (SGP4 + WGS-84)",
         "No refraction."],
        ["Off Boresight Angle η (°)", "angle between −r<sub>sat</sub> (nadir) and r<sub>GS</sub> − r<sub>sat</sub>",
         "The satellite antenna is assumed to point at nadir."],
        ["Doppler Shift (kHz)", M(r"f_D = -\dfrac{\dot{r}}{c}\, f, \quad \dot{r} = \mathbf{v}_{\mathrm{rel}} \cdot \hat{\mathbf{r}}"),
         "Downlink frequency; positive when the satellite approaches."],
        ["Doppler Rate (Hz/s)", [M(r"\dot{f}_D = -\dfrac{f\, \ddot{r}}{c}"),
                                 M(r"\ddot{r} = \dfrac{|\mathbf{v}_{\mathrm{rel}}|^2 - \dot{r}^2}{r} + \hat{\mathbf{r}} \cdot (\mathbf{a}_{\mathrm{sat}} - \mathbf{a}_{\mathrm{GS}})")], "a<sub>sat</sub> = −μ r<sub>sat</sub> / |r<sub>sat</sub>|³, "
                                                    "a<sub>GS</sub> = −ω<sub>E</sub>² (x, y, 0): exact at each "
                                                    "sample, also at the closest approach."],
        ["Sun Separation (°)", "angle between the Sun and the satellite seen from the station",
         M(r"\cos\theta = \sin e_1 \sin e_2 + \cos e_1 \cos e_2 \cos(a_1 - a_2)")],
    ], W3)
    story.append(P(
        "The <b>static</b> link budget has no orbit: it takes a circular orbit at the entered altitude h and a "
        "spherical Earth (<i>calculations.slant_range_from_elevation</i>, <i>off_nadir_angle_deg</i>):"))
    story.append(EQ(
        (r"d = \sqrt{(R + h)^2 - (R + h_{\mathrm{GS}})^2 \cos^2 \mathrm{el}} - (R + h_{\mathrm{GS}}) \sin \mathrm{el}",
         "Slant range."),
        (r"\eta = \arcsin\left[\dfrac{(R + h_{\mathrm{GS}}) \cos \mathrm{el}}{R + h}\right]",
         "Angle from the satellite nadir to the station."),
    ))

    # ----------------------------------------------------------------------
    story += H2("Atmospheric attenuation")
    story.append(P(
        "The slant-path attenuation of ITU-R P.618 is computed by the itur library "
        "(<i>itur.atmospheric_attenuation_slant_path</i>) at the station latitude, longitude and altitude, for "
        "the carrier frequency f, the elevation el, the percentage of time p and the antenna diameter D "
        "(scintillation aperture averaging, antenna efficiency 0.5):"))
    story.append(EQ(
        (r"p = 100 - \mathrm{availability}\ [\%]", "Static link budget: limited to 0.001–50 %."),
        (r"A = A_g + \sqrt{(A_r + A_c)^2 + A_s^2}", "ITU-R P.618-13 total slant-path attenuation."),
    ))
    story += BUL([
        "A<sub>g</sub>: gas, ITU-R P.676, with the water vapour density and content of the ITU-R maps at "
        "max(1, p) %.",
        "A<sub>c</sub>: cloud, ITU-R P.840 at max(1, p) %.",
        "A<sub>r</sub>: rain exceeded p % of the time, ITU-R P.618 with R<sub>0.01</sub> of P.837 and the rain "
        "height of P.839.",
        "A<sub>s</sub>: tropospheric scintillation, ITU-R P.618, only when Tropospheric scintillation is on.",
    ])
    story += BUL([
        "<b>Static link budget</b> (<i>atmospheric_attenuation_contributions</i>): <b>Clear Sky</b> column = "
        "A<sub>g</sub> only; <b>Rain Faded</b> column = A with gas, cloud, rain and, when switched on, "
        "scintillation, at the elevation of the tab.",
        "<b>Pass analysis</b> (<i>atmospheric_attenuation_vs_elevation</i>): A (rain faded, at p) is computed "
        "once on the elevation grid 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 75, 90° "
        "and interpolated linearly at the elevation of each sample (below 0.5° the 0.5° value). Every sample "
        "therefore carries the attenuation exceeded p % of the time: the pass results are rain-faded, not "
        "clear-sky.",
        "If the ITU-R model fails (e.g. outside its maps), the value at 5° is used at every elevation.",
        "The ITU-R tropospheric scintillation method is given for el ≥ 5° and about 4–20 GHz; outside it the "
        "results show a note (it is off by default).",
    ])

    # ----------------------------------------------------------------------
    story += H2("Ground station G/T")
    story.append(P(
        "Inputs (Downlink tab, shared with the Static LB · Downlink tab):"))
    story += TABLE([
        ["Parameter", "Key / unit", "Empty", "Use"],
        ["G/T GS", "gt_gs_dbk [dB/K]", "required", "G/T of the reference antenna at the reference elevation, "
                                                    "clear sky."],
        ["G/T Reference Diameter", "gt_ref_diameter_m [m]", "no scaling", "D<sub>ref</sub> of the antenna "
                                                                          "the G/T refers to."],
        ["G/T Reference Elevation", "gt_ref_elevation_deg [deg]", "5", "Elevation of G/T GS (and of Tsys)."],
        ["System Noise Temperature", "tsys_ref_k [K]", "constant G/T", "Tsys<sub>ref</sub> at the reference "
                                                                        "elevation, clear sky. Enables the "
                                                                        "elevation/weather model, the radome "
                                                                        "loss and the Sun noise."],
        ["Radome Loss", "radome_loss_db [dB]", "0", "Subtracted from the G/T (only with a System Noise "
                                                    "Temperature)."],
    ], W4, mono_cols=(1,))
    story.append(P("The G/T used at a sample or in a column is obtained in four steps "
                   "(<i>link_core.station_gt</i>, <i>calculations.ground_station_gt_dbk</i>):"))
    story.append(EQ(
        (r"(G/T)_{\mathrm{ref}} = (G/T)_{\mathrm{GS}} + 20\,\log_{10} \dfrac{D}{D_{\mathrm{ref}}}",
         "Step 1, without the station's own G/T (no scaling without D<sub>ref</sub>); otherwise "
         "(G/T)<sub>ref</sub> = gt_dbk of the catalogue."),
        (r"T_{\mathrm{sky}}(A) = T_m \left(1 - 10^{-A/10}\right) + T_c\, 10^{-A/10}",
         "Step 2: sky noise through an attenuation A [dB]; T<sub>m</sub> = 275 K, T<sub>c</sub> = 2.7 K."),
        (r"T_{\mathrm{sky,ref}} = T_{\mathrm{sky}}\!\left(A_g(\mathrm{el}_{\mathrm{ref}})\right)",
         "Clear sky (gas only, p = 50 %) at the reference elevation."),
        (r"T_{\mathrm{sky}} = T_{\mathrm{sky}}(A_g + A_c + A_r)",
         "At the sample or column; scintillation radiates no noise and is excluded."),
        (r"T_{\mathrm{sys}} = T_{\mathrm{sys,ref}} - T_{\mathrm{sky,ref}} + T_{\mathrm{sky}}",
         "Step 3: receiver and feed noise kept, sky noise replaced."),
        (r"G/T = (G/T)_{\mathrm{ref}} - 10\,\log_{10} \dfrac{T_{\mathrm{sys}}}{T_{\mathrm{sys,ref}}} - L_{\mathrm{radome}}", ""),
        (r"G/T \leftarrow G/T - \Delta_{\mathrm{Sun}}", "Step 4, pass analysis only (section §sun)."),
    ))
    story.append(P(
        "<b>Why the G/T changes with the elevation.</b> The G/T GS is a measured or specified value: the antenna "
        "gain G over the system noise temperature Tsys<sub>ref</sub> at one elevation, the reference elevation, "
        "in clear sky. The gain of the dish does not depend on where it points, but Tsys does: it is the noise "
        "of the receiver and feed, which is fixed, plus the noise of the sky the antenna looks at, T<sub>sky</sub>. "
        "An atmosphere that attenuates the signal by A dB also radiates noise (the absorbed power is re-emitted "
        "as heat at about Tm = 275 K): the longer the path through the atmosphere, the higher A and "
        "T<sub>sky</sub>. At low elevation the path is longest; rain and clouds add attenuation and noise too. "
        "The program therefore keeps the fixed part of Tsys, Tsys<sub>ref</sub> − T<sub>sky,ref</sub>, and adds "
        "the sky noise of the actual elevation and weather (step 3). The G/T changes by the ratio of the two "
        "system noise temperatures, as G is unchanged."))
    story += BUL([
        "At the reference elevation in clear sky T<sub>sky</sub> = T<sub>sky,ref</sub>, so G/T = G/T<sub>ref</sub> "
        "− radome loss: the radome loss is subtracted at every elevation, G/T GS being the value without radome.",
        "Below the reference elevation, or in rain, T<sub>sky</sub> &gt; T<sub>sky,ref</sub> and the G/T is lower; "
        "above it the G/T is slightly higher.",
        "The effect is larger for a low-noise receiver (small Tsys<sub>ref</sub>) and at higher frequencies.",
        "<b>The reference elevation matters</b>: the same G/T GS specified at zenith instead of 5° or 10° "
        "describes a worse receiver (the low zenith sky noise leaves more of Tsys<sub>ref</sub> to the "
        "receiver), so the G/T at low elevation comes out lower. When it is not known, 5° (the usual "
        "specification) is assumed; the station catalogue can give each station its own (gt_elevation_deg).",
        "Without a System Noise Temperature, steps 2–4 are skipped: G/T = G/T<sub>ref</sub> at every elevation "
        "and weather (no radome loss, no Sun noise).",
        "The static tab refuses a Tsys<sub>ref</sub> not higher than T<sub>sky,ref</sub>. Ground spill-over of the "
        "side lobes (ground noise picked up at low elevation) is not modelled.",
        "The static link budget computes it at its elevation, once for each column (clear: gas only; rain "
        "faded: gas + cloud + rain at p). The pass analysis computes it with the rain-faded attenuation on the "
        "elevation grid of section §atm and interpolates it at the elevation of each sample.",
    ])
    try:
        sky_ref, sky_ref_zenith, gt_rows = _gt_example()
        e = EXAMPLE
        story.append(NOTE(
            f"<b>Example</b> (invented station at {e['lat']:g}° N, {e['lon']:g}° E, {e['alt_m']:g} m, "
            f"{e['d_gs']:g} m dish, {e['freq_ghz']:g} GHz): G/T GS = {e['gt_dbk']:g} dB/K and Tsys<sub>ref</sub> = "
            f"{e['tsys_k']:g} K at {e['gt_ref_elevation_deg']:g}° clear sky, radome {e['radome_db']:g} dB, "
            f"availability {e['availability_pct']:g} % (p = {100 - e['availability_pct']:.1f} %). The clear-sky "
            f"gas attenuation at {e['gt_ref_elevation_deg']:g}° gives T<sub>sky,ref</sub> = {sky_ref:.1f} K, so "
            f"the receiver and feed account for {e['tsys_k'] - sky_ref:.1f} K. Computed by the program while this "
            "manual was built:"))
        table = [["El [°]", "A clear [dB]", "T_sky clear [K]", "Tsys clear [K]", "G/T clear [dB/K]",
                  "A rain [dB]", "T_sky rain [K]", "G/T rain [dB/K]", "G/T clear, ref. at zenith"]]
        for el, a_c, t_c, tsys_c, gt_c, a_r, t_r, gt_r, gt_z in gt_rows:
            table.append([f"{el:g}", f"{a_c:.3f}", f"{t_c:.1f}", f"{tsys_c:.1f}", f"{gt_c:.2f}", f"{a_r:.3f}",
                          f"{t_r:.1f}", f"{gt_r:.2f}", f"{gt_z:.2f}"])
        story += TABLE(table, [1.2, 1.8, 1.9, 1.9, 1.9, 1.7, 1.9, 1.9, 2.3], bold_cols=(0,))
        el5 = next(r for r in gt_rows if r[0] == 5.0)
        story.append(P(
            f"For example at 5° in clear sky: T<sub>sky</sub> = 275 (1 − 10<super>−{el5[1]:.3f}/10</super>) + 2.7 · "
            f"10<super>−{el5[1]:.3f}/10</super> = {el5[2]:.1f} K; Tsys = {e['tsys_k']:g} − {sky_ref:.1f} + "
            f"{el5[2]:.1f} = {el5[3]:.1f} K; G/T = {e['gt_dbk']:g} − 10·log({el5[3]:.1f} / {e['tsys_k']:g}) − "
            f"{e['radome_db']:g} = {el5[4]:.2f} dB/K. The last column takes the same {e['gt_dbk']:g} dB/K and "
            f"{e['tsys_k']:g} K as if specified at zenith (T<sub>sky,ref</sub> = {sky_ref_zenith:.1f} K): "
            f"{el5[8]:.2f} dB/K at 5°, {el5[4] - el5[8]:.2f} dB less."))
    except Exception as exc:  # pragma: no cover
        story.append(P(f"<i>[example not computed: {exc}]</i>"))

    # ----------------------------------------------------------------------
    story += H2("Sun noise")
    story.append(P(
        "Active in the pass analysis when <b>Sun noise</b> is ticked (key <font face='Courier'>sun_noise</font>, "
        "empty: on) and a System Noise Temperature is given (<i>calculations.sun_noise_temperature_k</i>). "
        "<b>Solar Flux</b> (key <font face='Courier'>solar_flux_sfu</font>, sfu = 10⁻²² W/m²/Hz, empty: quiet "
        "Sun) is the flux density at the downlink frequency."))
    story.append(EQ(
        (r"n = \mathrm{JD} - 2451545.0", "Days from J2000.0 (Astronomical Almanac, about 0.01°)."),
        (r"L = 280.460^{\circ} + 0.9856474^{\circ}\, n, \quad g = 357.528^{\circ} + 0.9856003^{\circ}\, n",
         "Mean longitude and mean anomaly of the Sun."),
        (r"\lambda_{\odot} = L + 1.915^{\circ} \sin g + 0.020^{\circ} \sin 2g, \quad \varepsilon = 23.439^{\circ} - 0.0000004^{\circ}\, n",
         "Ecliptic longitude and obliquity."),
        (r"\alpha = \mathrm{atan2}(\cos\varepsilon \sin\lambda_{\odot},\ \cos\lambda_{\odot}), \quad \delta = \arcsin(\sin\varepsilon \sin\lambda_{\odot})",
         ""),
        (r"\mathrm{GMST} = 280.46061837^{\circ} + 360.98564736629^{\circ}\, n, \quad H = \mathrm{GMST} + \lambda_{\mathrm{GS}} - \alpha",
         "Hour angle; elevation and azimuth from H, δ and the station latitude (no refraction)."),
        (r"T_b = \dfrac{S \cdot 10^{-22}\, \lambda^2}{2\, k\, \Omega_{\odot}}", "Solar Flux S [sfu] given."),
        (r"T_b = 120000\, f^{-0.75}\ \mathrm{K}", "Quiet Sun, f in GHz."),
        (r"\theta_{3\mathrm{dB}} = 70\, \lambda / D, \quad G_{\mathrm{peak}} = (G/T)_{\mathrm{ref}} + 10\,\log_{10} T_{\mathrm{sys,ref}}",
         "θ: Sun separation from the antenna axis."),
        (r"G(\theta) = G_{\mathrm{peak}} - 12 \left(\dfrac{\theta}{\theta_{3\mathrm{dB}}}\right)^2",
         "Main lobe, θ &lt; θ<sub>3dB</sub>."),
        (r"G(\theta) = \min\left[\max\left(G_{\mathrm{main}},\ G_{\mathrm{side}}\right),\ G_{\mathrm{peak}}\right], \quad G_{\mathrm{side}} = \max\left(32 - 25\,\log_{10} \theta,\ -10\right)",
         "θ ≥ θ<sub>3dB</sub>: ITU-R S.465 side lobes (θ taken ≥ 1°)."),
        (r"\Delta T = \dfrac{T_b\, \Omega_{\odot}\, G(\theta)}{4\pi}",
         "G linear; 0 with the Sun below −0.2665° of elevation."),
        (r"\Delta_{\mathrm{Sun}} = 10\,\log_{10}\left(1 + \dfrac{\Delta T}{T_{\mathrm{sys}}}\right)",
         "Sun Noise (dB): G/T and Eb/No drop by this value."),
    ))
    try:
        import calculations as c

        tb = c.sun_brightness_temperature_k(EXAMPLE["freq_ghz"])
        story.append(NOTE(
            f"At {EXAMPLE['freq_ghz']:g} GHz the quiet Sun has T<sub>b</sub> = {tb:,.0f} K. A pass is marked ☀ SUN NOISE in the "
            "contact list when the Sun noise reaches 0.5 dB; without a System Noise Temperature the list shows "
            "instead the passes with the Sun less than 5° from the satellite."))
    except Exception:  # pragma: no cover
        pass

    # ----------------------------------------------------------------------
    story += H2("Downlink parameters")
    story.append(P(
        "Fields of the Downlink tab (also shown, with blue labels, in the Static LB · Downlink tab). R<sub>c</sub> "
        "is the coded bit rate, OH the overhead, F the formatting factor of the static tab."))
    story += TABLE([
        ["Parameter", "Key / unit", "Empty", "Use"],
        ["Frequency", "frequency_ghz [GHz]", "required", "f: free space loss, ITU-R attenuation, Doppler, "
                                                         "λ for the beamwidth and the Sun noise, S4 scaling."],
        ["EIRP SAT", "eirp_sat_dbw [dBW]", "required", "EIRP at the peak of the satellite antenna pattern, used "
                                                       "when the static tab's EIRP source is 'EIRP value': "
                                                       "EIRP(η) = EIRP SAT − (G<sub>peak</sub> − G(η))."],
        ["Coded Bit Rate", "bitrate_mbps [Mbps]", "required", "R<sub>c</sub>, channel bit rate, coding "
                                                              "included."],
        ["Overhead (coded / info)", "overhead", "required", "OH = R<sub>c</sub> / R<sub>info</sub>: "
                                                            "R<sub>info</sub> = R<sub>c</sub> / OH."],
        ["Roll-off Factor", "rolloff", "no PFD check", "α: occupied bandwidth B = R<sub>c</sub> (1 + α) / "
                                                       "η<sub>s</sub>."],
        ["Spectral Efficiency", "spectral_efficiency_bpshz [b/s/Hz]", "no PFD check", "η<sub>s</sub> (QPSK: 2)."],
        ["Demodulator Loss", "demod_loss_db [dB]", "required", "Subtracted from C/No."],
        ["Link Availability", "link_availability_pct [%]", "required", "p = 100 − availability for rain and "
                                                                       "cloud."],
        ["Other Attenuations", "other_attenuations_db [dB]", "0", "Subtracted from the received power."],
        ["C/Io", "c_io_dbhz [dBHz]", "none", "Pass analysis only: C/(No+Io) = −10·log(10<super>−C/No/10"
                                                      "</super> + 10<super>−C/Io/10</super>)."],
        ["G/T GS and the G/T model", "see §gt", "", ""],
        ["Antenna Diameter GS", "antenna_diameter_m [m]", "3", "D of a station without diameter in the "
                                                               "catalogue."],
        ["Solar Flux, Sun noise", "solar_flux_sfu, sun_noise", "quiet Sun, on", "See §sun."],
        ["Tropospheric scintillation", "(switch)", "off", "Adds A<sub>s</sub> to the attenuation (§atm)."],
        ["Include the losses entered in the Static LB · Downlink tab", "(switch)", "on",
         "Takes from the static tab: the Tx chain, the ionospheric, polarisation and multipath losses, the "
         "modulation degradation, the formatting factor and the PFD limits (§static)."],
    ], W4, mono_cols=(1,))

    story += H2("Uplink parameters")
    story += TABLE([
        ["Parameter", "Key / unit", "Empty", "Use"],
        ["Uplink Frequency", "uplink_frequency_ghz [GHz]", "no uplink", "f<sub>UL</sub>: uplink free space loss "
                                                                        "and attenuation."],
        ["EIRP GS", "eirp_gs_dbw [dBW]", "no uplink", "Ground station EIRP (the station tracks: no pointing "
                                                      "loss)."],
        ["G/T SAT", "gt_sat_dbk [dB/K]", "no uplink", "Satellite G/T at the antenna pattern peak; the pattern "
                                                      "roll-off G<sub>peak</sub> − G(η) is subtracted."],
        ["Uplink Bit Rate", "uplink_bitrate_mbps [Mbps]", "downlink R<sub>c</sub>", "R<sub>c,UL</sub>."],
        ["Uplink Overhead", "uplink_overhead", "downlink OH", "OH<sub>UL</sub>."],
        ["Uplink Roll-off, Spectral Efficiency", "uplink_rolloff, uplink_spectral_efficiency_bpshz", "",
         "Only shown: information bit rate and bandwidth of the uplink in the Uplink tab."],
        ["Demodulator Loss SAT", "demod_loss_sat_db [dB]", "0", "Subtracted from the uplink C/No."],
        ["Uplink Other Attenuations", "uplink_other_attenuations_db [dB]", "0", "Subtracted from the uplink "
                                                                                 "received power."],
    ], W4, mono_cols=(1,))
    story.append(P("The uplink is computed only when Uplink Frequency, EIRP GS and G/T SAT are all given."))

    # ----------------------------------------------------------------------
    story += H2("Pass analysis: the columns of a sample")
    story.append(P(
        "For every sample above the horizon (<i>calculations.calculate_link_budget_parameters</i>; C/(No+Io) "
        "and Eb/No are meaningful in contact). L<sub>link</sub> = ionospheric + polarisation + multipath "
        "losses of the static tab, D<sub>mod</sub> its modulation degradation, F its formatting factor (each 0, "
        "0 and 1 when the static losses are not included). T<sub>chain</sub> = 10·log P<sub>Tx</sub> − circuit "
        "loss − VSWR loss of the static tab."))
    story += TABLE([
        ["Column", "Formula", "Notes"],
        ["Path Loss (dB)", M(r"L_{\mathrm{FS}} = 20\,\log_{10} d + 20\,\log_{10} f + 92.45"),
         "d in km, f in GHz; = 20 log(4π d / λ)."],
        ["Atmospheric Att (dB)", "A at the sample elevation", "§atm, rain faded at p."],
        ["S/C Antenna Gain (dBi)", "G(η): antenna_pattern.txt interpolated linearly", "0 dBi outside the table."],
        ["S/C Antenna Loss (dB)", M(r"\max\left[0,\ G_{\mathrm{peak}} - G(\eta)\right]"),
         "G<sub>peak</sub>: highest gain of the pattern."],
        ["EIRP (dBW)", [M(r"T_{\mathrm{chain}} + G(\eta)"),
                        M(r"\mathrm{EIRP}_{\mathrm{SAT}} - \left[G_{\mathrm{peak}} - G(\eta)\right]")],
         "First: static tab EIRP source 'Tx chain'; second: 'EIRP value'."],
        ["G/T GS (dB/K)", "§gt at the sample elevation, minus the Sun noise", ""],
        ["Link Losses (dB)", "L<sub>link</sub>", "No receive pointing loss: the station tracks."],
        ["Rx Power (dBW)", M(r"P_{\mathrm{Rx}} = \mathrm{EIRP} - L_{\mathrm{FS}} - A - L_{\mathrm{other}} - L_{\mathrm{link}}"), ""],
        ["C/(No+Io) (dBHz)", [M(r"C/N_0 = P_{\mathrm{Rx}} + G/T + 228.6 - L_{\mathrm{demod}} - D_{\mathrm{mod}}"),
                              M(r"C/(N_0 + I_0) = -10\,\log_{10}\left(10^{-\dfrac{C/N_0}{10}} + 10^{-\dfrac{C/I_0}{10}}\right)")],
         "The second line only when C/Io is given."],
        ["Eb/No (dB)", M(r"E_b/N_0 = C/(N_0 + I_0) - 10\,\log_{10} \dfrac{R_c\, F}{OH}"),
         "Referred to the transfer frame bit rate."],
        ["Doppler Shift, Doppler Rate", "§geom", "The contact list gives the largest |shift| and |rate| of each "
                                                "pass."],
        ["PFD (dBW/m²/4kHz)", M(r"\mathrm{EIRP} - 10\,\log_{10}(4\pi d^2) - 10\,\log_{10} \dfrac{B}{4\,\mathrm{kHz}}"), "d in m. When the static tab "
                                                                                     "has a PFD limit and α, "
                                                                                     "η<sub>s</sub> are set."],
        ["PFD Limit, PFD Margin", "limit(el) (§derived); margin = limit − PFD", "The contact list gives the "
                                                                             "smallest margin of each pass."],
        ["UL Path Loss, UL Atmospheric Att", "as above at f<sub>UL</sub>", ""],
        ["UL Rx Power (dBW)", [M(r"P_{\mathrm{Rx,UL}} = \mathrm{EIRP}_{\mathrm{GS}} - L_{\mathrm{FS,UL}} - A_{\mathrm{UL}} - L_{\mathrm{other,UL}}"),
                               M(r"\qquad - \left[G_{\mathrm{peak}} - G(\eta)\right] - L_{\mathrm{link,UL}}")],
         "Same satellite antenna pattern."],
        ["UL C/No (dBHz)", M(r"P_{\mathrm{Rx,UL}} + (G/T)_{\mathrm{SAT}} + 228.6 - L_{\mathrm{demod,SAT}} - D_{\mathrm{mod,UL}}"), ""],
        ["UL Eb/No (dB)", M(r"C/N_{0,\mathrm{UL}} - 10\,\log_{10} \dfrac{R_{c,\mathrm{UL}}\, F_{\mathrm{UL}}}{OH_{\mathrm{UL}}}"), ""],
    ], W3)

    # ----------------------------------------------------------------------
    story += H2("Static link budget inputs")
    story.append(P(
        "Keys of the sections <font face='Courier'>fixed_elevation_downlink</font> and "
        "<font face='Courier'>fixed_elevation_uplink</font> of parameters.json (written by <b>Save…</b>). The "
        "shared fields (frequency, EIRP, G/T, bit rate, overhead, roll-off, demodulator loss, other "
        "attenuations, availability, diameter) are the ones of §dl and §ul."))
    story += TABLE([
        ["Field", "Key / unit", "Empty", "Use"],
        ["Ground Station", "station", "", "Position for the ITU-R maps, the geomagnetic latitude and the G/T of "
                                         "the catalogue (independent of the pass analysis station)."],
        ["Elevation Angle", "elevation_deg [deg]", "5", "el of the budget (worst case: the minimum elevation)."],
        ["Satellite Altitude", "sat_altitude_km [km]", "required", "h: slant range and off-nadir angle."],
        ["Tropospheric scintillation", "scintillation", "off", "A<sub>s</sub> in the Rain Faded column."],
        ["Data framing", "framing_mode", "manual", "'ccsds' derives OH and F (§derived)."],
        ["RS Interleaving Depth I", "ccsds_interleaving", "5", "1–8."],
        ["Frame OCF / FECF / Packet Error Control", "ccsds_ocf / ccsds_fecf / ccsds_packet_error_control",
         "on / off / on", "4 B, 2 B, 2 B."],
        ["Source Packet Length", "ccsds_packet_length [B]", "not deducted", "L<sub>pkt</sub>."],
        ["Data Field Header (CUC)", "ccsds_packet_sec_header [B]", "0", "Packet secondary header."],
        ["Formatting Factor", "formatting_overhead", "1", "F = transfer frame bit rate / information bit rate."],
        ["Required Eb/No", "required_ebno_db [dB]", "no margin", "Threshold for the Data Recovery Margin."],
        ["Required Nominal Margin", "required_margin_db [dB]", "3", "Requirement on the Data Recovery Margin."],
        ["EIRP source", "eirp_mode", "EIRP value", "'tx_chain' or 'eirp'."],
        ["Transmitter Power", "tx_power_w [W]", "", "P<sub>Tx</sub> (Tx chain)."],
        ["Antenna Circuit Loss", "ant_circuit_loss_db [dB]", "0", "Cables, filters (Tx chain)."],
        ["VSWR", "vswr [:1]", "no loss", "Mismatch loss (Tx chain)."],
        ["Antenna Gain", "ant_gain_dbi [dBi]", "", "Tx antenna gain (Tx chain)."],
        ["Antenna gain from the pattern (DL)", "ant_gain_from_pattern", "off", "Antenna Gain = G(η) at this "
                                                                              "elevation; the Tx pointing loss "
                                                                              "is then 0."],
        ["Tx Antenna Axial Ratio", "ant_axial_ratio_db [dB]", "", "AR<sub>Tx</sub> for the polarisation loss."],
        ["Tx Antenna Pointing Loss", "tx_pointing_loss_db [dB]", "0", "Subtracted from Rx Power."],
        ["Ionospheric loss mode", "ionospheric_mode", "manual", "'s4' computes it from S4 (§derived)."],
        ["Ionospheric Loss", "ionospheric_loss_db [dB]", "0", "Manual value."],
        ["Ionospheric S4 (reference)", "iono_s4_ref", "0.3", "S4 at the reference frequency."],
        ["S4 Reference Frequency", "iono_s4_ref_freq_mhz [MHz]", "1575.42", "f<sub>ref</sub> of S4."],
        ["Polarisation Mismatch Loss", "polarisation_loss_db [dB]", "from the axial ratios, else 0",
         "Subtracted from Rx Power."],
        ["Multipath Losses", "multipath_loss_db [dB]", "0", "Subtracted from Rx Power."],
        ["PFD Limit 0–5° El", "pfd_limit [dBW/m²/4kHz]", "no PFD check", "Limit up to 5°."],
        ["PFD Limit 25–90° El", "pfd_limit_high [dBW/m²/4kHz]", "0–5° value", "Limit from 25°."],
        ["G/T SAT from the pattern (UL)", "gt_from_pattern", "off", "G/T SAT − (G<sub>peak</sub> − G(η)); the "
                                                                     "Rx depointing is then not used."],
        ["Modulation Degradation", "modulation_degradation_db [dB]", "0", "Added to the demodulator loss."],
        ["Rx Antenna Axial Ratio", "rx_axial_ratio_db [dB]", "", "AR<sub>Rx</sub>."],
        ["Rx Antenna Depointing Angle", "rx_depointing_deg [deg]", "no loss", "θ<sub>p</sub>."],
        ["Rx Antenna 3dB Beamwidth", "rx_beamwidth_deg [deg]", "70 λ / D (DL)", "θ<sub>3dB</sub> (required "
                                                                                "on the uplink with a "
                                                                                "depointing angle)."],
        ["Tolerances", "tolerances.&lt;key&gt;.favourable_db, adverse_db, distribution", "none",
         "Keys tx_power, tx_gain, propagation, rx_gt, rx_degradation (§derived)."],
    ], W4, mono_cols=(1,))

    # ----------------------------------------------------------------------
    story += H2("Static link budget: formulas of each result")
    story.append(P(
        "<i>calculations.calculate_fixed_elevation_link_budget</i> computes every row twice, Clear Sky and Rain "
        "Faded (§atm); C/Io is not used in the static budget."))
    story += TABLE([
        ["Row", "Formula", "Notes"],
        ["Tx Power [dBW]", M(r"10\,\log_{10} P_{\mathrm{Tx}}"), "P<sub>Tx</sub> in W; Tx chain."],
        ["VSWR Loss [dB]", M(r"-10\,\log_{10}(1 - \Gamma^2), \quad \Gamma = \dfrac{\mathrm{VSWR} - 1}{\mathrm{VSWR} + 1}"),
         "VSWR 1.5 → 0.177 dB."],
        ["Effective Gain [dBi]", M(r"G_{\mathrm{eff}} = G_{\mathrm{ant}} - L_{\mathrm{circuit}} - L_{\mathrm{VSWR}}"), ""],
        ["EIRP [dBW]", M(r"P_{\mathrm{Tx}} + G_{\mathrm{eff}}"), "Tx chain; otherwise the EIRP field."],
        ["Slant Range [km], Free Space Loss [dB]", M(r"L_{\mathrm{FS}} = 20\,\log_{10} d + 20\,\log_{10} f + 92.45"),
         "d from §geom, in km; f in GHz."],
        ["Gas / Cloud / Rain / Scintillation, Total Atmospheric Attenuation", "§atm", "Clear: gas only."],
        ["Ionospheric Loss [dB]", "manual, or S4 fade (§derived)", ""],
        ["Polarisation Mismatch Loss [dB]", "entered, or the co-polar value (§derived)", ""],
        ["Total Propagation Loss [dB]", M(r"L_{\mathrm{prop}} = L_{\mathrm{FS}} + A + L_{\mathrm{iono}} + L_{\mathrm{pol}} + L_{\mathrm{mp}} + L_{\mathrm{other}}"), ""],
        ["Rx Antenna Pointing Loss [dB]", M(r"12 \left(\dfrac{\theta_p}{\theta_{3\mathrm{dB}}}\right)^2"), "Parabolic main-lobe "
                                                                                    "approximation."],
        ["Rx Power [dBW]", M(r"P_{\mathrm{Rx}} = \mathrm{EIRP} - L_{\mathrm{prop}} - L_{\mathrm{point,Tx}} - L_{\mathrm{point,Rx}}"), ""],
        ["Rx G/T [dB/K]", "§gt for the condition (clear or rain sky noise)", "Uplink: G/T SAT (− roll-off)."],
        ["Received S/No [dBHz]", M(r"S/N_0 = P_{\mathrm{Rx}} + G/T + 228.6"), ""],
        ["C/No after Degradations [dBHz]", M(r"C/N_0 = S/N_0 - L_{\mathrm{demod}} - L_{\mathrm{mod}}"), ""],
        ["Transfer Frame Bit Rate [dBHz]", M(r"10\,\log_{10} R_{\mathrm{frame}}, \quad R_{\mathrm{frame}} = \dfrac{R_c\, F}{OH}"), ""],
        ["Received Eb/No [dB]", M(r"E_b/N_0 = C/N_0 - 10\,\log_{10} R_{\mathrm{frame}}"), ""],
        ["Data Recovery Margin [dB]", M(r"M = E_b/N_0 - (E_b/N_0)_{\mathrm{req}}"), "Green when ≥ Required Nominal "
                                                                         "Margin."],
        ["Margins with tolerances", "§derived", "MEAN − 3σ and Worst Case RSS ≥ 0 dB."],
        ["Power Flux Density [dBW/m²]", M(r"\mathrm{EIRP} - 10\,\log_{10}(4\pi d^2)"),
         "d in m; no atmospheric credit (conservative)."],
        ["Occupied Bandwidth [kHz]", M(r"B = \dfrac{R_c\,(1 + \alpha)}{\eta_s}"), ""],
        ["PFD in 4 kHz", M(r"\mathrm{PFD} - 10\,\log_{10} \dfrac{B}{4\,\mathrm{kHz}}"), "Uniform spectral density."],
        ["PFD Limit, PFD Margin", M(r"M_{\mathrm{PFD}} = \mathrm{limit}(\mathrm{el}) - \mathrm{PFD}_{4\mathrm{kHz}} \geq 0"), "limit(el): §derived."],
    ], W3)

    # ----------------------------------------------------------------------
    story += H2("Formulas of the derived losses and margins")
    story.append(H3("CCSDS TM framing (calculations.ccsds_tm_overheads)"))
    story.append(EQ(
        (r"L_{\mathrm{frame}} = 223\, I, \quad L_{\mathrm{cb}} = 255\, I", "Transfer frame and RS codeblock [bytes]."),
        (r"L_{\mathrm{df}} = L_{\mathrm{frame}} - 6 - L_{\mathrm{SH}} - 4_{\mathrm{OCF}} - 2_{\mathrm{FECF}}",
         "Data field: primary header 6 B, secondary header, OCF and FECF when selected."),
        (r"L_{\mathrm{user}} = L_{\mathrm{df}}\, \dfrac{L_{\mathrm{pkt}} - 6 - L_{\mathrm{CUC}} - L_{\mathrm{PEC}}}{L_{\mathrm{pkt}}}",
         "With a Source Packet length; otherwise L<sub>user</sub> = L<sub>df</sub>."),
        (r"F = \dfrac{L_{\mathrm{frame}}}{L_{\mathrm{user}}}, \quad c = 2\, \dfrac{L_{\mathrm{cb}} + 4}{L_{\mathrm{frame}}}",
         "Formatting and coding factors (ASM 4 B, convolutional 1/2)."),
        (r"OH = F\, c, \quad R_{\mathrm{frame}} = \dfrac{R_c}{c}, \quad R_{\mathrm{info}} = \dfrac{R_c}{OH}", ""),
    ))
    story.append(H3("Ionospheric scintillation (calculations.ionospheric_scintillation_fade_db)"))
    story.append(EQ(
        (r"S_4(f) = S_{4,\mathrm{ref}} \left(\dfrac{f}{f_{\mathrm{ref}}}\right)^{-1.5}", "ITU-R P.531, valid for S4 &lt; 0.4."),
        (r"P_{\mathrm{pp}} = 27.5\, S_4^{1.26}\ \mathrm{dB}, \quad L_{\mathrm{iono}} = P_{\mathrm{pp}} / 2",
         "Peak-to-peak fluctuation; the loss is half of it."),
        (r"\sin\phi_m = \sin\phi \sin\phi_p + \cos\phi \cos\phi_p \cos(\lambda - \lambda_p)",
         "Geomagnetic latitude; zones: equatorial |φ<sub>m</sub>| &lt; 20°, mid-latitude &lt; 58°, auroral "
         "oval &lt; 72°, polar cap above."),
    ))
    story.append(H3("Polarisation mismatch (calculations.polarisation_mismatch_loss_db)"))
    story.append(EQ(
        (r"r_1 = 10^{AR_{\mathrm{Tx}}/20}, \quad r_2 = 10^{AR_{\mathrm{Rx}}/20}", "Voltage axial ratios."),
        (r"L_{\mathrm{pol}} = -20\,\log_{10} \dfrac{2\,(r_1 + r_2)}{(1 + r_1)(1 + r_2)}", "Co-polar value (used)."),
        (r"L_{\mathrm{best,worst}} = -10\,\log_{10}\left[\dfrac{1}{2} + \dfrac{4 r_1 r_2 \pm (r_1^2 - 1)(r_2^2 - 1)}{2\,(1 + r_1^2)(1 + r_2^2)}\right]",
         "Range shown for information."),
    ))
    story.append(H3("PFD limit versus elevation (calculations.pfd_limit_at_elevation)"))
    story.append(EQ(
        (r"\mathrm{limit}(\mathrm{el}) = L_{\mathrm{low}}", "el ≤ 5°."),
        (r"\mathrm{limit}(\mathrm{el}) = L_{\mathrm{low}} + (L_{\mathrm{high}} - L_{\mathrm{low}})\, \dfrac{\mathrm{el} - 5^{\circ}}{20^{\circ}}",
         "5° &lt; el &lt; 25°."),
        (r"\mathrm{limit}(\mathrm{el}) = L_{\mathrm{high}}", "el ≥ 25°; without L<sub>high</sub>, L<sub>low</sub> at every elevation."),
    ))
    story.append(H3("Statistical margins (calculations.tolerance_statistics, margin_statistics)"))
    story.append(EQ(
        (r"a = -\mathrm{adv}, \quad b = +\mathrm{fav}", "Limits of each tolerance, in dB of margin."),
        (r"\mathrm{TRI}: \quad \mu = \dfrac{a + b}{3}, \quad \sigma^2 = \dfrac{a^2 + b^2 - ab}{18}", ""),
        (r"\mathrm{UNI}: \quad \mu = \dfrac{a + b}{2}, \quad \sigma^2 = \dfrac{(b - a)^2}{12}", ""),
        (r"\mathrm{GAU}: \quad \mu = \dfrac{a + b}{2}, \quad \sigma^2 = \left(\dfrac{b - a}{6}\right)^2", "Limits at ±3σ."),
        (r"M_{\mathrm{mean}} = M_{\mathrm{nom}} + \sum_k \mu_k, \quad \sigma = \sqrt{\sum_k \sigma_k^2}",
         "Independent contributions: means and variances add."),
        (r"M_{-3\sigma} = M_{\mathrm{mean}} - 3\sigma, \quad M_{\mathrm{RSS}} = M_{\mathrm{nom}} - \sqrt{\sum_k \mathrm{adv}_k^2}",
         "MEAN − 3σ and Worst Case RSS."),
        (r"M_{\mathrm{fav}} = M_{\mathrm{nom}} + \sum_k \mathrm{fav}_k, \quad M_{\mathrm{adv}} = M_{\mathrm{nom}} - \sum_k \mathrm{adv}_k",
         "All Favourable, All Adverse."),
    ))
    story.append(P(
        "Tolerances: Transmitter Power, Tx Antenna Gain, Atmospheric / Propagation Loss, Rx G/T, Receiver "
        "Degradation, entered as positive dB of effect on the margin; default distributions TRI, TRI, GAU, "
        "TRI, TRI. They are applied to the Rain Faded margin in the With Tolerances page."))

    # ----------------------------------------------------------------------
    story += H2("Worked example: static downlink")
    try:
        eta, gain, bandwidth, b, framing = _static_example()
        e = EXAMPLE
        clear, rain = b["clear"], b["rain_faded"]
        iono = b["ionospheric_scintillation"]
        pol = b["polarisation_loss_range_db"]
        tx = b["tx_power_dbw"]
        story.append(P(
            f"An invented S-band downlink: the station of §gt ({e['lat']:g}° N, {e['lon']:g}° E, {e['d_gs']:g} m), "
            f"el = {e['elevation_deg']:g}°, h = {e['altitude_km']:g} km, f = {e['freq_ghz']:g} GHz, availability "
            f"{e['availability_pct']:g} %; Tx chain {e['tx_power_w']:g} W, circuit loss {e['circuit_loss_db']:g} dB, "
            f"VSWR {e['vswr']:g}, gain from the satellite antenna pattern; coded rate {e['coded_mbps']:g} Mbps with "
            f"CCSDS framing I = {e['interleaving']}, OCF; S4 {e['s4_ref']:g} at GPS L1; axial ratios "
            f"{e['tx_axial_ratio_db']:g} / {e['rx_axial_ratio_db']:g} dB; multipath {e['multipath_db']:g} dB; Rx "
            f"depointing {e['rx_depointing_deg']:g}°; demodulator {e['demod_loss_db']:g} dB, modulation "
            f"{e['modulation_degradation_db']:g} dB; required Eb/No {e['required_ebno_db']:g} dB. Computed by the "
            "program while this manual was built; entered in the Static LB · Downlink tab they give the same "
            "values."))
        rows = [
            ["Step", "Calculation", "Value"],
            ["CCSDS framing", f"frame {framing['frame_bytes']} B, data field {framing['data_field_bytes']} B, "
                              f"codeblock {framing['codeblock_bytes']} B",
             f"F {framing['formatting']:.5f}, coding {framing['coding']:.5f}, OH {framing['total']:.5f}"],
            ["Off-nadir angle η", f"asin((R + h_GS) cos {e['elevation_deg']:g}° / (R + h))", f"{eta:.2f}°"],
            ["Antenna gain G(η)", "pattern interpolated", f"{gain:.2f} dBi"],
            ["Tx power", f"10·log {e['tx_power_w']:g} W", f"{tx:.2f} dBW"],
            ["VSWR loss", f"Γ = 0.3/2.3; −10·log(1 − Γ²)", f"{b['vswr_loss_db']:.3f} dB"],
            ["EIRP", f"{tx:.2f} − {e['circuit_loss_db']:g} − {b['vswr_loss_db']:.3f} "
                     f"{'+' if gain >= 0 else '−'} {abs(gain):.2f}",
             f"{b['eirp_dbw']:.2f} dBW"],
            ["Slant range d", "formula of §geom", f"{b['slant_range_km']:.1f} km"],
            ["Free space loss", f"20·log {b['slant_range_km']:.1f} + 20·log {e['freq_ghz']:g} + 92.45",
             f"{b['path_loss_db']:.2f} dB"],
            ["Atmosphere clear / rain", "gas / gas + cloud + rain",
             f"{clear['atmospheric_attenuation_db']:.3f} / {rain['atmospheric_attenuation_db']:.3f} dB"],
            ["Ionospheric loss", f"S4 {iono['s4_ref']:g} → {iono['s4']:.3f}; 27.5·S4^1.26 / 2",
             f"{b['ionospheric_loss_db']:.2f} dB"],
            ["Polarisation loss", f"co-polar, AR {e['tx_axial_ratio_db']:g} / {e['rx_axial_ratio_db']:g} dB",
             f"{pol['co_polar']:.3f} dB"],
            ["Rx pointing loss", f"12 ({e['rx_depointing_deg']:g} / {b['rx_beamwidth_3db_deg']:.2f})², "
                                 "θ3dB = 70 λ / D", f"{b['rx_pointing_loss_db']:.3f} dB"],
            ["Rx power clear", "EIRP − total propagation loss − pointing", f"{clear['rx_power_dbw']:.2f} dBW"],
            ["G/T clear", f"as in §gt at {e['elevation_deg']:g}°: Tsys {clear['system_noise_k']:.1f} K",
             f"{clear['gt_dbk']:.2f} dB/K"],
            ["S/No clear", "Rx power + G/T + 228.6", f"{clear['received_cno_dbhz']:.2f} dBHz"],
            ["C/No clear", f"− {e['demod_loss_db']:g} − {e['modulation_degradation_db']:g}",
             f"{clear['cno_dbhz']:.2f} dBHz"],
            ["Transfer frame bit rate", f"{e['coded_mbps']:g} Mbps / {framing['total']:.5f} × "
                                        f"{framing['formatting']:.5f}",
             f"{b['bitrate_formatted_bps'] / 1e3:.1f} kbps = {b['bitrate_dbhz']:.2f} dBHz"],
            ["Eb/No clear / rain", "C/No − 10·log R_frame", f"{clear['ebno_db']:.2f} / {rain['ebno_db']:.2f} dB"],
            ["Data Recovery Margin clear / rain", f"Eb/No − {e['required_ebno_db']:g}",
             f"{clear['margin_db']:.2f} / {rain['margin_db']:.2f} dB"],
            ["MEAN − 3σ (rain)", "Tx power ±0.5 dB TRI, propagation ±1 dB GAU",
             f"{rain['margin_statistics']['mean_minus_3sigma_db']:.2f} dB"],
            ["Worst Case RSS (rain)", "margin − √(0.5² + 1²)", f"{rain['margin_statistics']['worst_case_rss_db']:.2f} dB"],
            ["PFD in 4 kHz", f"B = {bandwidth / 1e6:g} MHz", f"{b['pfd_dbw_m2_per_4khz']:.2f} dBW/m²/4kHz"],
            ["PFD margin", f"limit {b['pfd_limit_dbw_m2_4khz']:g} at {e['elevation_deg']:g}°",
             f"{b['pfd_margin_db']:.2f} dB"],
        ]
        story += TABLE(rows, [4.0, 7.5, 5.0])
    except Exception as exc:  # pragma: no cover
        story.append(P(f"<i>[example not computed: {exc}]</i>"))

    # ======================================================================
    story += H1("Statistics reference")
    story.append(P(
        "The Statistics tab repeats the pass analysis for every day of the period, every station of the "
        "network and every satellite, then adds up the samples (<i>pass_statistics</i>). The per-sample "
        "values (elevation, contact, Eb/No) are the ones of chapter §ref, at 30 s steps."))

    story += H2("Inputs")
    story += TABLE([
        ["Field", "Key / unit", "Empty", "Use"],
        ["Start, Days", "[date], 1–366", "today, 7", "Period, from 00:00 UTC of the start date."],
        ["Ground stations", "statistics_stations [list]", "none (required)", "Network: a sample is usable "
                                                                             "through any of them."],
        ["Whole constellation", "(switch)", "on", "Custom-orbit mission: every satellite of §orbit; else the "
                                                  "satellite (plane, sat #) selected."],
        ["Required Eb/No", "(field) [dB]", "Static LB · Downlink value", "Eb/No<sub>req</sub>."],
        ["Required Margin", "(field) [dB]", "Static LB · Downlink value, else 0", "M<sub>req</sub>: a sample is "
                                                                                  "usable when Eb/No − "
                                                                                  "Eb/No<sub>req</sub> ≥ "
                                                                                  "M<sub>req</sub>."],
        ["Instrument Data, HKTM, Other Data", "downlink_data_kbps.instruments, .hktm, .other [kbps]",
         "not generated", "Data stored on board; their sum is the generation rate G."],
        ["Real-time TM", "downlink_data_kbps.realtime [kbps]", "0", "Sent live during the passes: the downlink "
                                                                   "rate for the stored data is R<sub>info</sub> − "
                                                                   "real-time."],
        ["On-board Memory", "memory_gb [GB]", "unlimited", "M<sub>max</sub> = GB × 8·10⁹ bit per satellite."],
        ["Min Pass Duration", "min_pass_min [min]", "0", "Shorter usable contacts are not used."],
        ["Antenna Setup", "antenna_setup_min [min]", "0", "Time to configure an antenna for another satellite."],
        ["Latency at end of pass", "latency_at_pass_end", "off", "Data count as on ground when the pass ends."],
        ["Latency requirement, Within", "latency_requirement.minutes, .percent", "none, 100 %", "Checked on "
                                                                                              "the downloaded "
                                                                                              "data."],
        ["Use backup antennas", "(switch)", "off", "Use all the antennas of each site (column antennas), not "
                                                   "one."],
    ], [2.8, 4.6, 2.9, 6.2], mono_cols=(1,))

    story += H2("Usable samples and passes")
    story.append(EQ(
        (r"m(s, g, t) = E_b/N_0\,(s, g, t) - (E_b/N_0)_{\mathrm{req}}", "s: satellite, g: station, t: sample."),
        (r"\mathrm{usable}(s, g, t) \Leftrightarrow \mathrm{contact}(s, g, t)\ \wedge\ m(s, g, t) \geq M_{\mathrm{req}}", ""),
        (r"C = R_d\, \Delta t, \quad R_d = \dfrac{R_c}{OH} - R_{\mathrm{RT}}",
         "Capacity of a usable sample; Δt = 30 s; R<sub>RT</sub>: real-time TM."),
    ))
    story.append(P(
        "<b>Min Pass Duration</b> (<i>drop_short_contacts</i>): each run of usable samples of one satellite over "
        "one station, followed across midnight, shorter than ⌈duration / 30 s⌉ samples is removed (set out of "
        "contact) before anything else is counted."))
    story.append(P(
        "<b>One satellite</b> (<i>combine_stations</i>, <i>aggregate</i>): at each sample the network is in "
        "contact when any station is; the station used is the one with the best margin. A pass is a run of "
        "samples in contact (one over midnight counts for the day it starts); its usable minutes are its "
        "usable samples × 0.5 min; the gap after a pass lasts until the next pass starts."))
    story.append(P(
        "<b>Several satellites</b> (<i>simulate_constellation</i>): every antenna serves one satellite at a "
        "time and every satellite uses one antenna. At each sample:"))
    story += BUL([
        "an antenna that served a satellite at the previous sample keeps it while that satellite stays usable "
        "there (a real antenna follows one pass);",
        "the free antennas are assigned by a maximum matching (augmenting paths), the satellites with more data "
        "on board first, each preferring its best-margin station; an antenna can take a satellite other than "
        "its last one only ⌈setup / 30 s⌉ samples after it stopped serving (Antenna Setup);",
        "a satellite that is usable somewhere but gets no antenna has a <b>conflict</b> for that sample;",
        "a site has 1 antenna, or its antennas column with Use backup antennas.",
    ])

    story += H2("On-board memory and latency")
    story.append(P(
        "Each satellite starts with an empty memory. At every sample t (<i>aggregate</i>, "
        "<i>simulate_constellation</i>):"))
    story.append(EQ(
        (r"Q \leftarrow Q + G\, \Delta t", "Q: data on board; G = instruments + HKTM + other."),
        (r"Q > Q_{\max}:\quad \mathrm{lost} \leftarrow \mathrm{lost} + Q - Q_{\max}, \quad Q \leftarrow Q_{\max}",
         "Memory full: the new data are lost."),
        (r"x = \min(Q,\ C)", "Bits downloaded in the sample, when the satellite is served (else 0)."),
        (r"Q \leftarrow Q - x", ""),
    ))
    story.append(P(
        "The data leave the memory first in, first out. With D the bits downloaded before this sample, the "
        "x bits downloaded now were generated between D / G and (D + x) / G seconds from the start, uniformly. They "
        "reach the ground at t<sub>arr</sub> = end of the sample, or, with <b>Latency at end of pass</b>, at the "
        "end of the pass (the run of served samples) they belong to. Each sample adds a chunk:"))
    story.append(EQ(
        (r"\ell_{\min} = t_{\mathrm{arr}} - \dfrac{D + x}{G}, \quad \ell_{\max} = t_{\mathrm{arr}} - \dfrac{D}{G}",
         "Latency range of the chunk of x bits."),
        (r"s(L) = \dfrac{\sum_k x_k\, \mathrm{clip}\left(\frac{60\,L - \ell_{\min,k}}{\ell_{\max,k} - \ell_{\min,k}},\ 0,\ 1\right)}{\sum_k x_k}",
         "Share of the data on ground within L minutes."),
        (r"L_X = \min\,\{L : s(L) \geq X\}", "Latency at X % (bisection, 60 steps)."),
        (r"\bar{\ell} = \dfrac{\sum_k x_k\, (\ell_{\min,k} + \ell_{\max,k})/2}{\sum_k x_k}, \quad \ell_{\mathrm{max}} = \max_k\, \ell_{\max,k}",
         "Mean and maximum latency."),
    ))
    story.append(P(
        "So <i>'80 % within 67 min'</i> means: of all the bits generated on board and downloaded in the period "
        "(all satellites together for a constellation), 80 % reached the ground at most 67 min after they were "
        "generated. It is a share of the data volume, not of the time or of the passes. Data still on board at "
        "the end of the period are not counted. When the memory overflows the latency of the rest is "
        "approximate."))

    story += H2("Results")
    story += TABLE([
        ["Result", "Definition", "Notes"],
        ["Passes, per day", "number of passes (service windows per satellite for a constellation) / days", ""],
        ["Contact time [min/day]", "samples in contact × 0.5 min / days", "Served samples for a constellation."],
        ["Max / mean gap [min]", "between the end of a pass and the start of the next", "One satellite."],
        ["Usable time [min/day]", "usable samples × 0.5 min / days", ""],
        ["Downlink capacity [MB/day]", "usable samples × R<sub>d</sub> · 30 s / 8·10⁶ per day (mean, min, max)",
         ""],
        ["Generated [MB/day]", "G · 86400 / 8·10⁶ (× satellites)", ""],
        ["Downloaded [MB/day]", "Σ x / 8·10⁶ per day", "≤ capacity, ≤ generated."],
        ["On board [MB]", "memory at the end of each sample (sum of all satellites); max and final", ""],
        ["Memory use [%]", "100 · largest memory of one satellite / M<sub>max</sub>", ""],
        ["Autonomy [h]", "M<sub>max</sub> / G / 3600: time to fill the memory without contacts", ""],
        ["Lost [MB, %]", "Σ lost / 8·10⁶; 100 · lost / (G · period)", ""],
        ["Latency mean / max / at X %, within L", "§lat above", "Requirement met when share(L) ≥ X %."],
        ["Per satellite: served, conflicts [min/day]", "served / conflict samples × 0.5 min / days", ""],
        ["Per station: served [min/day], volume [MB/day]", "samples served by its antennas; data downloaded "
                                                          "through it", ""],
        ["Per station: busy [%]", "100 · served samples / (samples × antennas)", ""],
        ["Per station: contention [min/day]", "samples with more usable satellites than antennas × 0.5 min / "
                                              "days", ""],
    ], [4.5, 8.5, 3.5])

    story += H2("Worked example: memory and latency")
    try:
        runs = _statistics_example()
        s, s_end = runs[False].summary, runs[True].summary
        p50, p80 = runs[False].latency_percentile(50), runs[False].latency_percentile(80)
        story.append(P(
            "One satellite generating G = 60 kbps, information rate R<sub>d</sub> = 1.5 Mbps (45 Mbit per "
            "sample), two contacts of 10 min from 60 to 70 min and from 160 to 170 min, 4 hours. Computed with "
            "<i>pass_statistics.aggregate</i> while this manual was built:"))
        story += BUL([
            "At 60 min the memory holds 60 kbps × 3630 s = 217.8 Mbit (27.2 MB, the first sample's generation "
            "included); 45 Mbit leave per sample, so the backlog is gone in 5 samples (2.5 min), then each "
            "sample sends only the 1.8 Mbit generated in it.",
            "The oldest bit, generated at t = 0, lands at the end of the first contact sample (60.5 min): "
            "latency 60.5 min. The data generated during the contact land within a sample.",
            f"Result: mean latency {s['mean_latency_min']:.1f} min, max {s['max_latency_min']:.1f} min, "
            f"50 % within {p50:.1f} min, 80 % within {p80:.1f} min; {s['latency_within_pct']:.1f} % within "
            "75 min.",
            f"With Latency at end of pass every chunk lands at 70 min (first pass) or 170 min (second): mean "
            f"{s_end['mean_latency_min']:.1f} min, max {s_end['max_latency_min']:.1f} min.",
        ])
    except Exception as exc:  # pragma: no cover
        story.append(P(f"<i>[example not computed: {exc}]</i>"))
    return story


def parameters_json_rows():
    """Rows (key, label, unit, description) of every parameters.json key of the schema."""
    import mission_schema

    rows = []
    group = None
    for g, path, label, unit, _kind, description in mission_schema.SCHEMA:
        rows.append([g if g != group else "", ".".join(path), label + (f" [{unit}]" if unit else ""), description])
        group = g
    return rows
