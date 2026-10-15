"""Chapters of the user manual that document every parameter and how it is computed.

Called by build_manual.py, which passes its own module (``m``) for the
building blocks (H1, H2, P, TABLE, ...) so the chapter and figure numbers
stay in one sequence. The worked examples are computed here, at build time,
with the functions of the program, so the numbers printed in the manual are
the ones the program gives.
"""

import datetime
import json
import os
import warnings

import numpy as np


def _fmt(value, digits=2):
    return f"{value:.{digits}f}"


# --------------------------------------------------------------------------
# Worked examples (computed with the program's own functions)
# --------------------------------------------------------------------------
def _static_example(root):
    """AWS downlink over Svalbard at 3°, as the Static LB · Downlink tab computes it."""
    import astropy.units as u

    import calculations as c

    with open(os.path.join(root, "missions", "AWS", "parameters.json"), encoding="utf-8") as f:
        p = json.load(f)
    dl = p["fixed_elevation_downlink"]
    c.load_antenna_pattern(os.path.join(root, "missions", "AWS", "antenna_pattern.txt"))
    lat, lon, alt_m = c.GROUND_STATIONS["Svalbard"]
    eta = c.off_nadir_angle_deg(dl["elevation_deg"], dl["sat_altitude_km"], alt_m / 1000.0)
    gain = c.satellite_antenna_gain_dbi(eta)
    rate = p["bitrate_mbps"] * 1e6
    bandwidth = rate * (1 + p["rolloff"]) / p["spectral_efficiency_bpshz"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        budget = c.calculate_fixed_elevation_link_budget(
            freq=p["frequency_ghz"] * u.GHz, elevation_deg=dl["elevation_deg"],
            sat_altitude_km=dl["sat_altitude_km"], lat_gs=lat, lon_gs=lon, alt_gs_km=alt_m / 1000.0,
            d_gs=c.GROUND_STATION_DIAMETERS.get("Svalbard", p["antenna_diameter_m"]), eirp=p["eirp_sat_dbw"],
            gt=p["gt_gs_dbk"], demod_loss=p["demod_loss_db"], bitrate=rate, overhead=p["overhead"],
            other_att=p["other_attenuations_db"] or 0.0, pointing_loss_db=0.0,
            link_availability_pct=p["link_availability_pct"], include_scintillation=bool(dl.get("scintillation")),
            required_ebno=dl["required_ebno_db"], tx_power_w=dl["tx_power_w"],
            antenna_circuit_loss_db=dl["ant_circuit_loss_db"], vswr=dl["vswr"], antenna_gain_dbi=gain,
            ionospheric_s4_ref=dl["iono_s4_ref"], ionospheric_s4_ref_freq_mhz=dl["iono_s4_ref_freq_mhz"],
            polarisation_loss_db=dl["polarisation_loss_db"], multipath_loss_db=dl["multipath_loss_db"],
            modulation_degradation_db=dl["modulation_degradation_db"], occupied_bandwidth_hz=bandwidth,
            pfd_limit_dbw_m2_4khz=dl["pfd_limit"], pfd_limit_high_dbw_m2_4khz=dl["pfd_limit_high"],
            tx_axial_ratio_db=dl["ant_axial_ratio_db"], rx_axial_ratio_db=dl["rx_axial_ratio_db"],
            rx_depointing_deg=dl["rx_depointing_deg"], rx_beamwidth_3db_deg=dl["rx_beamwidth_deg"],
            formatting_overhead=dl["formatting_overhead"],
            tolerances=[dict(name="Tx Antenna Gain", key=k, **v) for k, v in dl["tolerances"].items()],
            gt_model={"ref_elevation_deg": p["gt_ref_elevation_deg"], "tsys_ref_k": p["tsys_ref_k"],
                      "radome_loss_db": p["radome_loss_db"]},
        )
    framing = c.ccsds_tm_overheads(int(dl["ccsds_interleaving"]), dl["ccsds_ocf"], dl["ccsds_fecf"])
    return p, dl, eta, gain, bandwidth, budget, framing


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
def chapters(m, root):
    """The reference chapters, as a list of flowables."""
    # Section references are written as tokens (§gt, §atm, ...) and numbered here, so they follow the
    # position of these chapters in the manual.
    ref_ch = m._chapter[0] + 1
    sections = ["const", "orbit", "stations", "geom", "atm", "gt", "sun", "dl", "ul", "columns", "static",
                "static_results", "derived", "example"]
    tokens = {name: f"{ref_ch}.{k}" for k, name in enumerate(sections, start=1)}
    tokens.update(ref=str(ref_ch), lat=f"{ref_ch + 1}.3", calc=str(m.CALC_CHAPTER[0]))

    def T(text):
        for name in sorted(tokens, key=len, reverse=True):
            text = text.replace("§" + name, tokens[name])
        return text

    H1, H2, H3, CODE = m.H1, m.H2, m.H3, m.CODE

    def P(text):
        return m.P(T(text))

    def NOTE(text):
        return m.NOTE(T(text))

    def BUL(items):
        return m.BUL([T(item) for item in items])

    def TABLE(rows, *args, **kwargs):
        return m.TABLE([[T(str(cell)) for cell in row] for row in rows], *args, **kwargs)

    def EQ(*lines):
        text = T("\n".join(lines)).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        return CODE(text.replace("  ", "&nbsp; "))

    W4 = [3.0, 4.6, 1.7, 7.2]  # parameter tables: name, key / unit, empty value, use
    W3 = [4.2, 4.8, 7.5]  # result tables: name, formula, notes

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
        "n  = √(μ / a³)                      p = a (1 − e²)",
        "SSO inclination:  cos i = − (dΩ/dt)_SSO / (1.5 · n · J2 · (R / p)²)",
        "RAAN = α_sun(epoch) + 15° · (LTAN − 12 h)          (α_sun: right ascension of the Sun)",
        "mean anomaly M = u − ω              (u: argument of latitude of the satellite at the epoch)",
        "k = 0.75 · J2 · (R / p)²",
        "dM/dt = n · (1 + k · √(1 − e²) · (3 cos² i − 1))      dω/dt = n · k · (5 cos² i − 1)",
        "nodal period T = 2π / (dM/dt + dω/dt)",
    ))
    story.append(P(
        "The epoch is 00:00 UTC of the day analysed (of the first day for the statistics); the drag term is 0. "
        "The mean motion written in the TLE is adjusted (up to 5 iterations, to 10⁻¹⁰) until SGP4's own "
        "dM/dt + dω/dt equals 2π / T: SGP4 then flies the J2 nodal period, and the pass times stay right over the "
        "analysed days (a Keplerian mean motion would drift by about a minute per orbit). The Sun right "
        "ascension uses the low-precision formulae of the Astronomical Almanac (section §sun)."))
    story.append(P(
        "<b>Constellation</b> (<i>link_core.constellation_members</i>): with P planes in the LTAN list, the "
        "satellite j = 1…N of plane p = 0…P−1 has the argument of latitude at the epoch"))
    story.append(EQ("u = (j − 1) · 360° / N + p · phasing    (mod 360°)"))
    story.append(P(
        "An extra satellite \"LTAN@u\" is added with that plane and that u, numbered after the regular "
        "satellites of its plane. The satellites are named \"&lt;mission&gt; &lt;LTAN&gt; #j\" and get the "
        "catalogue numbers 99900, 99901, …"))
    try:
        import tle_source as ts

        alt, ecc = 595.6, 0.001214
        inc = ts.sso_inclination_deg(alt, ecc)
        period = ts.nodal_period_min(ts.EARTH_RADIUS_KM + alt, ecc, 97.794)
        story.append(NOTE(
            f"<b>Example – EPS-Sterna</b> (595.6 km, e = 0.001214): the sun-synchronous inclination is "
            f"{inc:.3f}° (the parameters give the reference value 97.794°, which is used); the nodal period at "
            f"97.794° is {period:.2f} min. With planes 15:30, 19:30, 23:30, N = 2 and phasing 60°, u is 0° and "
            f"180° (15:30), 60° and 240° (19:30), 120° and 300° (23:30): the six satellites are 60° apart."))
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
        "required(t) = max(Min Elevation, mask(azimuth(t)))",
        "in contact  ⇔  elevation(t) ≥ required(t)",
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
        ["Doppler Shift (kHz)", "f<sub>D</sub> = −ṙ / c · f, ṙ = v<sub>rel</sub> · r̂",
         "Downlink frequency; positive when the satellite approaches."],
        ["Doppler Rate (Hz/s)", "−f · r̈ / c, r̈ = (|v<sub>rel</sub>|² − ṙ²) / r + r̂ · (a<sub>sat</sub> − "
                                "a<sub>GS</sub>)", "a<sub>sat</sub> = −μ r<sub>sat</sub> / |r<sub>sat</sub>|³, "
                                                    "a<sub>GS</sub> = −ω<sub>E</sub>² (x, y, 0): exact at each "
                                                    "sample, also at the closest approach."],
        ["Sun Separation (°)", "angle between the Sun and the satellite seen from the station",
         "cos θ = sin e₁ sin e₂ + cos e₁ cos e₂ cos(a₁ − a₂)."],
    ], W3)
    story.append(P(
        "The <b>static</b> link budget has no orbit: it takes a circular orbit at the entered altitude h and a "
        "spherical Earth (<i>calculations.slant_range_from_elevation</i>, <i>off_nadir_angle_deg</i>):"))
    story.append(EQ(
        "d = √((R + h)² − ((R + h_GS) · cos el)²) − (R + h_GS) · sin el",
        "η = asin((R + h_GS) · cos el / (R + h))",
    ))

    # ----------------------------------------------------------------------
    story += H2("Atmospheric attenuation")
    story.append(P(
        "The slant-path attenuation of ITU-R P.618 is computed by the itur library "
        "(<i>itur.atmospheric_attenuation_slant_path</i>) at the station latitude, longitude and altitude, for "
        "the carrier frequency f, the elevation el, the percentage of time p and the antenna diameter D "
        "(scintillation aperture averaging, antenna efficiency 0.5):"))
    story.append(EQ(
        "p = 100 − Link Availability [%]          (static: limited to 0.001 … 50 %)",
        "A_g : gas, ITU-R P.676, water vapour density and content of the ITU-R maps at max(1, p) %",
        "A_c : cloud, ITU-R P.840 at max(1, p) %",
        "A_r : rain, ITU-R P.618 with R0.01 of P.837 and the rain height of P.839, exceeded p %",
        "A_s : tropospheric scintillation, ITU-R P.618 (only when Tropospheric scintillation is on)",
        "A   = A_g + √((A_r + A_c)² + A_s²)",
    ))
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
        "1. G/T_ref = gt_dbk of the station                         when the catalogue gives it",
        "           = G/T GS + 20·log(D / D_ref)                     otherwise (0 without D_ref)",
        "2. T_sky(A) = Tm · (1 − 10^(−A/10)) + Tc · 10^(−A/10)",
        "   T_sky,ref = T_sky(A_g at the reference elevation, p = 50 %, gas only)",
        "   T_sky     = T_sky(A_g + A_c + A_r of the sample or column)    (no scintillation: it radiates",
        "                                                               no noise)",
        "3. Tsys = Tsys_ref − T_sky,ref + T_sky",
        "   G/T  = G/T_ref − 10·log(Tsys / Tsys_ref) − radome loss",
        "4. pass analysis only: G/T −= Sun Noise (dB)                  (section §sun)",
    ))
    story.append(P(
        "Without a System Noise Temperature, steps 2–4 are skipped: G/T = G/T<sub>ref</sub> at every elevation "
        "and weather (no radome loss, no Sun noise). The static tab refuses a Tsys<sub>ref</sub> not higher than "
        "T<sub>sky,ref</sub>. Ground spill-over of the side lobes is not modelled. In the pass analysis the G/T is "
        "computed on the elevation grid of section §atm and interpolated."))

    # ----------------------------------------------------------------------
    story += H2("Sun noise")
    story.append(P(
        "Active in the pass analysis when <b>Sun noise</b> is ticked (key <font face='Courier'>sun_noise</font>, "
        "empty: on) and a System Noise Temperature is given (<i>calculations.sun_noise_temperature_k</i>). "
        "<b>Solar Flux</b> (key <font face='Courier'>solar_flux_sfu</font>, sfu = 10⁻²² W/m²/Hz, empty: quiet "
        "Sun) is the flux density at the downlink frequency."))
    story.append(EQ(
        "Sun position (Astronomical Almanac, ≈ 0.01°):  n = JD − 2451545.0",
        "  L = 280.460° + 0.9856474° n     g = 357.528° + 0.9856003° n",
        "  λ = L + 1.915° sin g + 0.020° sin 2g     ε = 23.439° − 0.0000004° n",
        "  α = atan2(cos ε sin λ, cos λ)    δ = asin(sin ε sin λ)",
        "  GMST = 280.46061837° + 360.98564736629° n     H = GMST + longitude − α",
        "  elevation and azimuth from H, δ and the station latitude (no refraction)",
        "Brightness temperature:  T_b = S · 10^-22 · λ² / (2 k Ω_sun)      (Solar Flux S given)",
        "                         T_b = 120000 · f[GHz]^-0.75 K           (quiet Sun)",
        "Antenna gain θ from the axis (θ = Sun separation):  θ3dB = 70 λ / D",
        "  G_peak = G/T_ref + 10·log(Tsys_ref)",
        "  θ < θ3dB : G(θ) = G_peak − 12 (θ / θ3dB)²",
        "  θ ≥ θ3dB : G(θ) = min(max(main lobe, side lobe), G_peak)",
        "  side lobe (ITU-R S.465) = max(32 − 25·log θ, −10) dBi,  θ taken ≥ 1°",
        "Added noise:  ΔT = T_b · Ω_sun · G(θ)_lin / 4π      (0 when the Sun elevation ≤ −0.2665°)",
        "Sun Noise (dB) = 10·log(1 + ΔT / Tsys)        G/T and Eb/No drop by this value",
    ))
    try:
        import calculations as c

        tb = c.sun_brightness_temperature_k(1.707)
        story.append(NOTE(
            f"At 1.707 GHz the quiet Sun has T<sub>b</sub> = {tb:,.0f} K. A pass is marked ☀ SUN NOISE in the "
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
        ["Path Loss (dB)", "FSL = 20·log d[km] + 20·log f[GHz] + 92.45", "= 20·log(4π d / λ)."],
        ["Atmospheric Att (dB)", "A at the sample elevation", "§atm, rain faded at p."],
        ["S/C Antenna Gain (dBi)", "G(η): antenna_pattern.txt interpolated linearly", "0 dBi outside the table."],
        ["S/C Antenna Loss (dB)", "G<sub>peak</sub> − G(η) (≥ 0)", "G<sub>peak</sub>: highest gain of the "
                                                                   "pattern."],
        ["EIRP (dBW)", "T<sub>chain</sub> + G(η)   (static tab EIRP source 'Tx chain')<br/>"
                       "EIRP SAT − (G<sub>peak</sub> − G(η))   ('EIRP value')", ""],
        ["G/T GS (dB/K)", "§gt at the sample elevation, minus the Sun noise", ""],
        ["Link Losses (dB)", "L<sub>link</sub>", "No receive pointing loss: the station tracks."],
        ["Rx Power (dBW)", "EIRP − FSL − A − Other Att. − L<sub>link</sub>", ""],
        ["C/(No+Io) (dBHz)", "C/No = Rx Power + G/T + 228.6 − Demod Loss − D<sub>mod</sub>; combined with "
                             "C/Io when given", ""],
        ["Eb/No (dB)", "C/(No+Io) − 10·log(R<sub>c</sub> / OH × F)", "Referred to the transfer frame bit rate."],
        ["Doppler Shift, Doppler Rate", "§geom", "The contact list gives the largest |shift| and |rate| of each "
                                                "pass."],
        ["PFD (dBW/m²/4kHz)", "EIRP − 10·log(4π) − 20·log d[m] − 10·log(B / 4000 Hz)", "When the static tab "
                                                                                     "has a PFD limit and α, "
                                                                                     "η<sub>s</sub> are set."],
        ["PFD Limit, PFD Margin", "limit(el) (§derived); margin = limit − PFD", "The contact list gives the "
                                                                             "smallest margin of each pass."],
        ["UL Path Loss, UL Atmospheric Att", "as above at f<sub>UL</sub>", ""],
        ["UL Rx Power (dBW)", "EIRP GS − FSL<sub>UL</sub> − A<sub>UL</sub> − UL Other Att. − (G<sub>peak</sub> − "
                              "G(η)) − L<sub>link,UL</sub>", "Same satellite antenna pattern."],
        ["UL C/No (dBHz)", "UL Rx Power + G/T SAT + 228.6 − Demod Loss SAT − D<sub>mod,UL</sub>", ""],
        ["UL Eb/No (dB)", "UL C/No − 10·log(R<sub>c,UL</sub> / OH<sub>UL</sub> × F<sub>UL</sub>)", ""],
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
        ["Tx Power [dBW]", "10·log P<sub>Tx</sub>[W]", "Tx chain."],
        ["VSWR Loss [dB]", "−10·log(1 − Γ²), Γ = (VSWR − 1) / (VSWR + 1)", "VSWR 1.5 → 0.177 dB."],
        ["Effective Gain [dBi]", "Antenna Gain − Circuit Loss − VSWR Loss", ""],
        ["EIRP [dBW]", "Tx Power + Effective Gain   (Tx chain), else the EIRP field", ""],
        ["Slant Range [km], Free Space Loss [dB]", "§geom; 20·log d + 20·log f + 92.45", ""],
        ["Gas / Cloud / Rain / Scintillation, Total Atmospheric Attenuation", "§atm", "Clear: gas only."],
        ["Ionospheric Loss [dB]", "manual, or S4 fade (§derived)", ""],
        ["Polarisation Mismatch Loss [dB]", "entered, or the co-polar value (§derived)", ""],
        ["Total Propagation Loss [dB]", "FSL + A + ionospheric + polarisation + multipath + other", ""],
        ["Rx Antenna Pointing Loss [dB]", "12 · (θ<sub>p</sub> / θ<sub>3dB</sub>)²", "Parabolic main-lobe "
                                                                                    "approximation."],
        ["Rx Power [dBW]", "EIRP − Total Propagation Loss − Tx pointing − Rx pointing", ""],
        ["Rx G/T [dB/K]", "§gt for the condition (clear or rain sky noise)", "Uplink: G/T SAT (− roll-off)."],
        ["Received S/No [dBHz]", "Rx Power + G/T + 228.6", ""],
        ["C/No after Degradations [dBHz]", "S/No − Receiver Degradation − Modulation Degradation", ""],
        ["Transfer Frame Bit Rate [dBHz]", "10·log(R<sub>c</sub> / OH × F)", ""],
        ["Received Eb/No [dB]", "C/No after Degradations − Transfer Frame Bit Rate", ""],
        ["Data Recovery Margin [dB]", "Received Eb/No − Required Eb/No", "Green when ≥ Required Nominal "
                                                                         "Margin."],
        ["Margins with tolerances", "§derived", "MEAN − 3σ and Worst Case RSS ≥ 0 dB."],
        ["Power Flux Density [dBW/m²]", "EIRP − 10·log(4π d²), d in m", "No atmospheric credit (conservative)."],
        ["Occupied Bandwidth [kHz]", "B = R<sub>c</sub> (1 + α) / η<sub>s</sub>", ""],
        ["PFD in 4 kHz", "PFD − 10·log(B / 4 kHz)", "Uniform spectral density."],
        ["PFD Limit, PFD Margin", "limit(el); limit − PFD in 4 kHz ≥ 0", "§derived."],
    ], W3)

    # ----------------------------------------------------------------------
    story += H2("Formulas of the derived losses and margins")
    story.append(H3("CCSDS TM framing (calculations.ccsds_tm_overheads)"))
    story.append(EQ(
        "frame = 223 · I bytes                                  codeblock = 255 · I bytes",
        "data field = frame − 6 (primary header) − secondary header − 4 (OCF) − 2 (FECF)",
        "user data  = data field                                       (no packet length)",
        "           = data field · (L_pkt − 6 − CUC − PEC) / L_pkt     (with a Source Packet length)",
        "F (formatting) = frame / user data",
        "coding = (codeblock + 4 (ASM)) / frame · 2 (convolutional 1/2)",
        "OH = F · coding        R_frame = R_c / coding        R_info = R_c / OH",
    ))
    story.append(H3("Ionospheric scintillation (calculations.ionospheric_scintillation_fade_db)"))
    story.append(EQ(
        "S4(f) = S4_ref · (f / f_ref)^-1.5                (ITU-R P.531, valid for S4 < 0.4)",
        "peak-to-peak fluctuation = 27.5 · S4^1.26 dB      ionospheric loss = half of it",
        "geomagnetic latitude: sin φm = sin φ sin φp + cos φ cos φp cos(λ − λp)",
        "zone: equatorial |φm| < 20°, mid-latitude < 58°, auroral oval < 72°, polar cap above",
    ))
    story.append(H3("Polarisation mismatch (calculations.polarisation_mismatch_loss_db)"))
    story.append(EQ(
        "r1 = 10^(AR_Tx/20)     r2 = 10^(AR_Rx/20)",
        "co-polar loss = −20·log(2 (r1 + r2) / ((1 + r1)(1 + r2)))      (used)",
        "best / worst = −10·log(0.5 + (4 r1 r2 ± (r1² − 1)(r2² − 1)) / (2 (1 + r1²)(1 + r2²)))  (shown)",
    ))
    story.append(H3("PFD limit versus elevation (calculations.pfd_limit_at_elevation)"))
    story.append(EQ(
        "limit(el) = L_low                                   el ≤ 5°",
        "          = L_low + (L_high − L_low) (el − 5) / 20   5° < el < 25°",
        "          = L_high                                  el ≥ 25°",
        "(without L_high: L_low at every elevation)",
    ))
    story.append(H3("Statistical margins (calculations.tolerance_statistics, margin_statistics)"))
    story.append(EQ(
        "for each tolerance: a = −adverse, b = +favourable (dB of margin)",
        "  TRI: mean (a + b) / 3   variance (a² + b² − ab) / 18",
        "  UNI: mean (a + b) / 2   variance (b − a)² / 12",
        "  GAU: mean (a + b) / 2   variance ((b − a) / 6)²       (limits at ±3σ)",
        "Margin Mean = nominal + Σ mean      σ = √(Σ variance)",
        "MEAN − 3σ = Margin Mean − 3σ        Worst Case RSS = nominal − √(Σ adverse²)",
        "All Favourable = nominal + Σ favourable     All Adverse = nominal − Σ adverse",
    ))
    story.append(P(
        "Tolerances: Transmitter Power, Tx Antenna Gain, Atmospheric / Propagation Loss, Rx G/T, Receiver "
        "Degradation, entered as positive dB of effect on the margin; default distributions TRI, TRI, GAU, "
        "TRI, TRI. They are applied to the Rain Faded margin in the With Tolerances page."))

    # ----------------------------------------------------------------------
    story += H2("Worked example: static downlink")
    try:
        p, dl, eta, gain, bandwidth, b, framing = _static_example(root)
        clear, rain = b["clear"], b["rain_faded"]
        iono = b["ionospheric_scintillation"]
        story.append(P(
            "AWS, Static LB · Downlink as saved in missions/AWS/parameters.json: Svalbard, el = "
            f"{dl['elevation_deg']:g}°, h = {dl['sat_altitude_km']:g} km, f = {p['frequency_ghz']:g} GHz, "
            f"availability {p['link_availability_pct']:g} %. Computed by the program while this manual was "
            "built; the Static LB tab shows the same values (rounded)."))
        tx = b["tx_power_dbw"]
        rows = [
            ["Step", "Calculation", "Value"],
            ["Off-nadir angle η", f"asin((R + h_GS) cos {dl['elevation_deg']:g}° / (R + h))", f"{eta:.2f}°"],
            ["Antenna gain G(η)", "antenna_pattern.txt", f"{gain:.2f} dBi"],
            ["Tx power", f"10·log {dl['tx_power_w']:g} W", f"{tx:.2f} dBW"],
            ["VSWR loss", f"VSWR {dl['vswr']:g}", f"{b['vswr_loss_db']:.3f} dB"],
            ["EIRP", f"{tx:.2f} − {dl['ant_circuit_loss_db']:g} − {b['vswr_loss_db']:.3f} + {gain:.2f}",
             f"{b['eirp_dbw']:.2f} dBW"],
            ["Slant range d", "formula of §geom", f"{b['slant_range_km']:.1f} km"],
            ["Free space loss", f"20·log {b['slant_range_km']:.1f} + 20·log {p['frequency_ghz']:g} + 92.45",
             f"{b['path_loss_db']:.2f} dB"],
            ["Atmosphere clear / rain", "gas / gas + cloud + rain",
             f"{clear['atmospheric_attenuation_db']:.3f} / {rain['atmospheric_attenuation_db']:.3f} dB"],
            ["Ionospheric loss", f"S4 {iono['s4_ref']:g} → {iono['s4']:.3f}; 27.5·S4^1.26 / 2",
             f"{b['ionospheric_loss_db']:.2f} dB"],
            ["Polarisation loss", "entered", f"{b['polarisation_loss_db']:.2f} dB"],
            ["Rx pointing loss", f"12 ({dl['rx_depointing_deg']:g} / {dl['rx_beamwidth_deg']:g})²",
             f"{b['rx_pointing_loss_db']:.3f} dB"],
            ["Rx power clear", "EIRP − total propagation loss − pointing", f"{clear['rx_power_dbw']:.2f} dBW"],
            ["Sky noise ref. / clear", f"gas at {p['gt_ref_elevation_deg']:g}° / at {dl['elevation_deg']:g}°",
             f"{b['gt_model']['sky_ref_k']:.1f} / {clear['sky_noise_k']:.1f} K"],
            ["Tsys clear", f"{p['tsys_ref_k']:g} − {b['gt_model']['sky_ref_k']:.1f} + {clear['sky_noise_k']:.1f}",
             f"{clear['system_noise_k']:.1f} K"],
            ["G/T clear", f"{p['gt_gs_dbk']:g} − 10·log({clear['system_noise_k']:.1f}/{p['tsys_ref_k']:g}) − "
                          f"{p['radome_loss_db']:g}", f"{clear['gt_dbk']:.2f} dB/K"],
            ["S/No clear", "Rx power + G/T + 228.6", f"{clear['received_cno_dbhz']:.2f} dBHz"],
            ["C/No clear", f"− {p['demod_loss_db']:g} − {dl['modulation_degradation_db']:g}",
             f"{clear['cno_dbhz']:.2f} dBHz"],
            ["Transfer frame bit rate", f"{p['bitrate_mbps']:g} Mbps / {p['overhead']:g} × "
                                        f"{dl['formatting_overhead']:g}",
             f"{b['bitrate_formatted_bps'] / 1e3:.1f} kbps = {b['bitrate_dbhz']:.2f} dBHz"],
            ["Eb/No clear / rain", "C/No − 10·log R_frame", f"{clear['ebno_db']:.2f} / {rain['ebno_db']:.2f} dB"],
            ["Data Recovery Margin clear / rain", f"Eb/No − {dl['required_ebno_db']:g}",
             f"{clear['margin_db']:.2f} / {rain['margin_db']:.2f} dB"],
            ["MEAN − 3σ (rain)", "Tx gain ±1.5 dB TRI: σ = 1.5/√6",
             f"{rain['margin_statistics']['mean_minus_3sigma_db']:.2f} dB"],
            ["PFD in 4 kHz", f"B = {bandwidth / 1e6:g} MHz", f"{b['pfd_dbw_m2_per_4khz']:.2f} dBW/m²/4kHz"],
            ["PFD margin", f"limit {b['pfd_limit_dbw_m2_4khz']:g} at {dl['elevation_deg']:g}°",
             f"{b['pfd_margin_db']:.2f} dB"],
            ["CCSDS I = 5", "frame 1115 B, data field 1105 B, codeblock 1275 B",
             f"F {framing['formatting']:.5f}, coding {framing['coding']:.5f}, OH {framing['total']:.5f}"],
        ]
        story += TABLE(rows, [4.0, 7.5, 5.0])
        story.append(P(
            "The same link was checked against SatMaster Pro (validation/satmaster_aws_svalbard.md): C/No in "
            "clear sky within 0.15 dB."))
    except Exception as exc:  # pragma: no cover - only if the mission files are missing
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
        "s: satellite, g: station, t: sample",
        "margin(s, g, t) = Eb/No(s, g, t) − Eb/No_req",
        "usable(s, g, t) ⇔ in contact(s, g, t) and margin ≥ M_req",
        "capacity of a usable sample = R_d · 30 s,   R_d = R_c / OH − real-time [bit/s]",
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
        "memory += G · 30 s                                        (G = instruments + HKTM + other)",
        "if memory > M_max: lost += memory − M_max, memory = M_max",
        "take = min(memory, capacity of the sample)   if the satellite is served, else 0",
        "memory −= take",
    ))
    story.append(P(
        "The data leave the memory first in, first out. With D the bits downloaded before this sample, the "
        "bits taken now were generated between D / G and (D + take) / G seconds from the start, uniformly. They "
        "reach the ground at t<sub>arr</sub> = end of the sample, or, with <b>Latency at end of pass</b>, at the "
        "end of the pass (the run of served samples) they belong to. Each sample adds a chunk:"))
    story.append(EQ(
        "chunk = (bits = take,  latency from t_arr − (D + take)/G  to  t_arr − D/G)",
        "share within L minutes = Σ over chunks of bits · clip((60 L − lat_min) / (lat_max − lat_min), 0, 1) / Σ bits",
        "latency at X % = smallest L with share(L) ≥ X %      (bisection, 60 steps)",
        "mean latency = Σ bits · (lat_min + lat_max) / 2 / Σ bits        max latency = largest lat_max",
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
        ["Downloaded [MB/day]", "Σ take / 8·10⁶ per day", "≤ capacity, ≤ generated."],
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
