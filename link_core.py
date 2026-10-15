"""Pass analysis without the user interface.

The window (gui.py), the multi-day statistics and the command line
(lbtool.py) all compute the passes here, with the same models:

* :class:`Station` and :class:`LinkInputs` hold the inputs;
* :func:`static_link_values` turns a static link budget section of
  ``parameters.json`` (or the values of a Static LB tab) into the losses the
  pass analysis takes from it, and :func:`pass_loss_kwargs` into keyword
  arguments of :func:`calculations.calculate_link_budget_parameters`;
* :class:`PassAnalysis` computes one day of samples (geometry, visibility
  with the minimum elevation and horizon mask, link budget with the
  atmospheric attenuation and G/T at each elevation, PFD) and the contact
  windows.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Any, Optional

import astropy.units as u
import numpy as np
import pandas as pd
from skyfield.api import EarthSatellite, load, wgs84

import calculations

STEP_S = 30


@dataclass
class Station:
    name: str
    lat: float
    lon: float
    alt_m: float
    d_gs: float = 3.0
    # Minimum elevation (deg) per azimuth, as (azimuths, elevations), or None.
    mask: Optional[tuple] = None

    @property
    def alt_km(self) -> float:
        return self.alt_m / 1000.0


@dataclass
class UplinkInputs:
    freq_ghz: float
    eirp_gs: float
    gt_sat: float
    demod_loss: Optional[float] = None
    other_att: Optional[float] = None
    bitrate_bps: Optional[float] = None
    overhead: Optional[float] = None


@dataclass
class LinkInputs:
    """Downlink (and optional uplink) inputs of the pass analysis."""

    freq_ghz: float
    link_availability_pct: float
    eirp_sat: float
    gt_gs: float  # of the reference antenna (gt_ref_diameter_m), at the reference elevation
    demod_loss: float
    bitrate_bps: float  # coded
    overhead: float
    other_att: float = 0.0
    cio_db: Optional[float] = None
    scintillation: bool = False
    gt_ref_diameter_m: Optional[float] = None
    # {"ref_elevation_deg", "tsys_ref_k", "radome_loss_db"} or None for a constant G/T.
    gt_model: Optional[dict] = None
    loss_kwargs: dict = field(default_factory=lambda: {"apply_antenna_pattern": True})
    uplink: Optional[UplinkInputs] = None
    # PFD check: limits (dBW/m² in 4 kHz) and occupied bandwidth (Hz).
    pfd_limit_low: Optional[float] = None
    pfd_limit_high: Optional[float] = None
    occupied_bandwidth_hz: Optional[float] = None
    # Sun in the ground antenna beam: noise added to the system noise temperature
    # (needs gt_model). Solar flux at the downlink frequency in sfu; None: quiet Sun.
    sun_noise: bool = True
    solar_flux_sfu: Optional[float] = None

    @property
    def p(self) -> float:
        return 100.0 - self.link_availability_pct

    @property
    def info_bitrate_bps(self) -> float:
        return self.bitrate_bps / self.overhead

    def gt_for(self, station: Station) -> float:
        """G/T at the reference elevation of the station's antenna (diameter scaling)."""
        return self.gt_gs + calculations.gt_diameter_scaling_db(station.d_gs, self.gt_ref_diameter_m)


def _value(section: dict, key: str, default=None):
    value = section.get(key, default)
    return default if value in (None, "") else value


def static_link_values(section: dict, freq_ghz: float, d_gs: float) -> dict:
    """What the pass analysis takes from a static link budget section.

    ``section`` has the keys of ``fixed_elevation_downlink`` /
    ``fixed_elevation_uplink`` in parameters.json (numbers, or None when
    empty). Returns the Tx chain (dBW, before the antenna gain) when the
    EIRP comes from it, the elevation-independent losses, the modulation
    degradation, the formatting factor, the PFD limits and the required
    Eb/No and margin. Raises ValueError for an invalid value.
    """
    if section.get("ionospheric_mode") == "s4":
        iono = dict(
            ionospheric_s4_ref=float(section["iono_s4_ref"]),
            ionospheric_s4_ref_freq_mhz=float(section["iono_s4_ref_freq_mhz"]),
        )
    else:
        iono = dict(ionospheric_loss_db=float(_value(section, "ionospheric_loss_db", 0.0)))
    losses = calculations.static_link_losses(
        freq_ghz,
        d_gs,
        polarisation_loss_db=_value(section, "polarisation_loss_db"),
        tx_axial_ratio_db=_value(section, "ant_axial_ratio_db"),
        rx_axial_ratio_db=_value(section, "rx_axial_ratio_db"),
        multipath_loss_db=float(_value(section, "multipath_loss_db", 0.0)),
        rx_depointing_deg=_value(section, "rx_depointing_deg"),
        rx_beamwidth_3db_deg=_value(section, "rx_beamwidth_deg"),
        **iono,
    )
    formatting = float(_value(section, "formatting_overhead", 1.0))
    if formatting <= 0:
        raise ValueError("'Formatting Factor (frame / user data)' must be greater than 0.")
    tx_chain_dbw = None
    if section.get("eirp_mode") == "tx_chain":
        power_w = _value(section, "tx_power_w")
        if power_w is None or float(power_w) <= 0:
            raise ValueError("'Transmitter Power [W]' must be greater than 0.")
        vswr = _value(section, "vswr")
        tx_chain_dbw = (
            10 * math.log10(float(power_w))
            - float(_value(section, "ant_circuit_loss_db", 0.0))
            - (calculations.vswr_mismatch_loss_db(float(vswr)) if vswr else 0.0)
        )
    return {
        "tx_chain_dbw": tx_chain_dbw,
        "ionospheric_loss_db": losses["ionospheric_loss_db"],
        "polarisation_loss_db": losses["polarisation_loss_db"],
        "multipath_loss_db": losses["multipath_loss_db"],
        "rx_pointing_loss_db": losses["rx_pointing_loss_db"],
        "tx_pointing_loss_db": float(_value(section, "tx_pointing_loss_db", 0.0)),
        "modulation_degradation_db": float(_value(section, "modulation_degradation_db", 0.0)),
        "formatting_overhead": formatting,
        "pfd_limit": _value(section, "pfd_limit"),
        "pfd_limit_high": _value(section, "pfd_limit_high"),
        "required_ebno_db": _value(section, "required_ebno_db"),
        "required_margin_db": _value(section, "required_margin_db"),
    }


def pass_loss_kwargs(downlink: Optional[dict], uplink: Optional[dict] = None) -> dict:
    """Keyword arguments of calculate_link_budget_parameters from :func:`static_link_values`.

    The satellite antenna pattern is always applied; the downlink EIRP is
    the Tx chain plus the pattern gain when the static tab uses the Tx
    chain. The ground station tracks the satellite: no pointing loss.
    ``None`` for a direction leaves out its static losses.
    """
    kwargs: dict = {"apply_antenna_pattern": True}
    for values, prefix in ((downlink, ""), (uplink, "ul_")):
        if values is None:
            continue
        if prefix == "" and values["tx_chain_dbw"] is not None:
            kwargs["tx_chain_dbw"] = values["tx_chain_dbw"]
        kwargs[f"{prefix}link_losses_db"] = (
            values["ionospheric_loss_db"] + values["polarisation_loss_db"] + values["multipath_loss_db"]
        )
        kwargs[f"{prefix}modulation_degradation_db"] = values["modulation_degradation_db"]
        kwargs[f"{prefix}formatting_overhead"] = values["formatting_overhead"]
    return kwargs


def pass_direction(sat, ts, start_time: datetime, end_time: datetime) -> str:
    """``A`` for an ascending pass, ``D`` for a descending one."""
    lat_start = wgs84.subpoint_of(sat.at(ts.from_datetime(start_time))).latitude.degrees
    lat_end = wgs84.subpoint_of(sat.at(ts.from_datetime(end_time))).latitude.degrees
    return "A" if lat_end >= lat_start else "D"


# Columns rounded in the samples table (decimals; 2 when not listed).
_ROUNDED = (
    "Elevation (°)", "Slant Range (km)", "Path Loss (dB)", "Atmospheric Att (dB)", "G/T GS (dB/K)",
    "S/C Antenna Gain (dBi)", "EIRP (dBW)", "S/C Antenna Loss (dB)", "Off Boresight Angle (°)",
    "Link Losses (dB)", "UL Link Losses (dB)", "Rx Power (dBW)", "C/(No+Io) (dBHz)", "Eb/No (dB)",
    "Doppler Shift (kHz)", "UL Path Loss (dB)", "UL Atmospheric Att (dB)", "UL S/C Antenna Loss (dB)",
    "UL Rx Power (dBW)", "UL C/No (dBHz)", "UL Eb/No (dB)",
)
_DECIMALS = {"Atmospheric Att (dB)": 3, "UL Atmospheric Att (dB)": 3}


class PassAnalysis:
    """One satellite seen from one station: samples of a day and contact windows.

    ``link=None`` gives the geometry and the contacts only (a satellite
    without link budget parameters).
    """

    def __init__(
        self,
        tle: tuple[str, str],
        station: Station,
        link: Optional[LinkInputs] = None,
        min_elevation_deg: float = calculations.MIN_ELEVATION_DEG,
        step_s: int = STEP_S,
    ):
        self.ts = load.timescale()
        self.sat = EarthSatellite(tle[0], tle[1], "SAT", self.ts)
        self.station = station
        self.gs = wgs84.latlon(station.lat, station.lon, elevation_m=station.alt_m)
        self.link = link
        self.min_elevation_deg = min_elevation_deg
        self.step_s = step_s
        self.freq = (link.freq_ghz if link else 1.0) * u.GHz
        self._tables: dict[str, Any] = {}

    # --- elevation tables (ITU-R model evaluated once, then interpolated) ---

    def _downlink_tables(self):
        if "dl" not in self._tables:
            link, st = self.link, self.station
            grid, attenuation = calculations.atmospheric_attenuation_vs_elevation(
                st.lat, st.lon, link.freq_ghz, link.p, st.d_gs, st.alt_km,
                include_scintillation=link.scintillation,
            )
            gt_ref = link.gt_for(st)
            tsys_grid = None
            if link.gt_model is None:
                gt_grid = np.full_like(attenuation, gt_ref)
            else:
                model = link.gt_model
                if model["tsys_ref_k"] <= 0:
                    raise ValueError("The system noise temperature must be greater than 0 K.")
                sky_ref = calculations.gt_reference_sky_noise_k(
                    st.lat, st.lon, link.freq_ghz, st.d_gs, st.alt_km, model["ref_elevation_deg"]
                )
                sky = calculations.sky_noise_temperature_k(attenuation)
                gt_grid = calculations.ground_station_gt_dbk(
                    gt_ref, model["tsys_ref_k"], sky_ref, sky, model.get("radome_loss_db", 0.0),
                )
                tsys_grid = model["tsys_ref_k"] - sky_ref + sky
            self._tables["dl"] = (grid, attenuation, gt_grid, tsys_grid)
        return self._tables["dl"]

    def _uplink_table(self):
        if "ul" not in self._tables:
            ul, st = self.link.uplink, self.station
            self._tables["ul"] = calculations.atmospheric_attenuation_vs_elevation(
                st.lat, st.lon, ul.freq_ghz, self.link.p, st.d_gs, st.alt_km,
                include_scintillation=self.link.scintillation,
            )
        return self._tables["ul"]

    def atmosphere_summary(self) -> str:
        """One line on the attenuation and G/T used, for the window."""
        grid, attenuation, gt_grid, _tsys = self._downlink_tables()
        low = calculations.MIN_ELEVATION_DEG
        text = (
            f"Atmospheric Att (dB): {calculations.interpolate_vs_elevation(low, grid, attenuation):.3f} "
            f"at {low:g}° El, {attenuation[-1]:.3f} at 90° (per sample)\n"
        )
        if self.link.gt_model is None:
            text += f"G/T GS {gt_grid[0]:g} dB/K at all elevations"
        else:
            text += (
                f"G/T GS {calculations.interpolate_vs_elevation(low, grid, gt_grid):.2f} dB/K at {low:g}°, "
                f"{gt_grid[-1]:.2f} at 90°"
            )
        scaling = calculations.gt_diameter_scaling_db(self.station.d_gs, self.link.gt_ref_diameter_m)
        if scaling:
            text += f" ({self.station.d_gs:g} m antenna: {scaling:+.2f} dB from the reference antenna)"
        if self._sun_applies():
            flux = self.link.solar_flux_sfu
            text += (
                f"\nSun noise: {'quiet Sun' if flux is None else f'{flux:g} sfu'}, "
                f"+{self._sun_noise_k(np.array([0.0]), np.array([90.0]))[0]:.0f} K on the antenna axis"
            )
        return text

    def _sun_applies(self) -> bool:
        return bool(self.link is not None and self.link.sun_noise and self.link.gt_model is not None)

    def _sun_noise_k(self, separation_deg, sun_elevation_deg):
        """Noise temperature (K) the Sun adds to the station antenna."""
        link = self.link
        peak_gain = link.gt_for(self.station) + 10 * np.log10(link.gt_model["tsys_ref_k"])
        return calculations.sun_noise_temperature_k(
            separation_deg, sun_elevation_deg, peak_gain, self.station.d_gs, link.freq_ghz, link.solar_flux_sfu
        )

    # --- one day ---

    def times(self, day_start: datetime) -> list[datetime]:
        return [day_start + timedelta(seconds=self.step_s * i) for i in range(86400 // self.step_s)]

    def day(self, day_start: datetime, budget: str = "above_horizon", times: Optional[list] = None) -> pd.DataFrame:
        """Samples of the day starting at ``day_start`` (UTC midnight).

        ``budget`` chooses where the link budget is computed: "above_horizon"
        (every sample above 0°, enough for the pass plots and tables),
        "visible" (only the samples in contact) or "all".
        """
        times = times or self.times(day_start)
        sky_times = self.ts.from_datetimes(times)
        altitudes, slant_ranges, dopplers, off_boresight = calculations.prepare_topocentric_data(
            self.sat, self.gs, sky_times, self.freq
        )
        azimuth = calculations.topocentric_azimuth_deg(self.sat, self.gs, sky_times)
        required = np.maximum(
            self.min_elevation_deg, calculations.mask_elevation_deg(azimuth, self.station.mask)
        )
        visible = np.asarray(altitudes) >= required
        sun_el, sun_az = calculations.sun_altaz_deg(self.station.lat, self.station.lon, times)
        sun_sep = calculations.angular_separation_deg(altitudes, azimuth, sun_el, sun_az)
        frame = pd.DataFrame({"Time (UTC)": times})
        if self.link is None:
            frame["Elevation (°)"] = np.round(altitudes, 2)
            frame["Slant Range (km)"] = np.round(slant_ranges, 2)
        else:
            frame = self._budget(frame, times, sky_times, altitudes, slant_ranges, dopplers, off_boresight,
                                 visible, budget, sun_sep, sun_el)
        frame["Sun Separation (°)"] = np.round(sun_sep, 2)
        frame["Azimuth (°)"] = np.round(azimuth, 2)
        frame["Min Elevation (°)"] = np.round(required, 2)
        frame["Visible"] = visible
        if self.link is not None:
            self._add_pfd(frame)
        return frame

    def _budget(self, frame, times, sky_times, altitudes, slant_ranges, dopplers, off_boresight, visible, budget,
                sun_sep, sun_el):
        link, st = self.link, self.station
        grid, attenuation, gt_grid, tsys_grid = self._downlink_tables()
        atm = calculations.interpolate_vs_elevation(np.asarray(altitudes, dtype=float), grid, attenuation)
        gt = calculations.interpolate_vs_elevation(np.asarray(altitudes, dtype=float), grid, gt_grid)
        sun_db = None
        if self._sun_applies():
            # The Sun raises the system noise temperature: G/T and Eb/No drop by the same dB.
            tsys = calculations.interpolate_vs_elevation(np.asarray(altitudes, dtype=float), grid, tsys_grid)
            sun_db = 10 * np.log10(1.0 + self._sun_noise_k(sun_sep, sun_el) / tsys)
            gt = gt - sun_db
        ul_atm = None
        ul = link.uplink
        if ul is not None:
            ul_grid, ul_attenuation = self._uplink_table()
            ul_atm = calculations.interpolate_vs_elevation(np.asarray(altitudes, dtype=float), ul_grid, ul_attenuation)
        if budget == "all":
            indices = range(len(times))
        elif budget == "visible":
            indices = np.flatnonzero(visible)
        else:
            indices = np.flatnonzero(np.asarray(altitudes) > 0.0)
        cisat_lin = 10 ** (link.cio_db / 10.0) if link.cio_db is not None else None
        rows = []
        for i in indices:
            params = calculations.calculate_link_budget_parameters(
                sky_times[i], self.sat, self.gs, self.freq, link.p, st.d_gs, st.alt_km,
                link.eirp_sat, float(gt[i]), link.demod_loss, link.bitrate_bps, link.overhead, cisat_lin,
                link.other_att, float(atm[i]), altitudes[i], slant_ranges[i], dopplers[i], off_boresight[i],
                ul.freq_ghz * u.GHz if ul else None,
                ul.eirp_gs if ul else None,
                ul.gt_sat if ul else None,
                ul.demod_loss if ul else None,
                ul.other_att if ul else None,
                float(ul_atm[i]) if ul_atm is not None else None,
                ul.bitrate_bps if ul else None,
                ul.overhead if ul else None,
                **link.loss_kwargs,
            )
            params["G/T GS (dB/K)"] = float(gt[i])
            if sun_db is not None:
                params["Sun Noise (dB)"] = round(float(sun_db[i]), 2)
            for key in _ROUNDED:
                if params.get(key) is not None:
                    params[key] = round(params[key], _DECIMALS.get(key, 2))
            params["Time (UTC)"] = times[i]
            rows.append(params)
        computed = pd.DataFrame(rows)
        if computed.empty:
            frame["Elevation (°)"] = np.round(altitudes, 2)
            frame["Slant Range (km)"] = np.round(slant_ranges, 2)
            return frame
        merged = frame.merge(computed, on="Time (UTC)", how="left")
        # Geometry everywhere, also where the budget was not computed.
        merged["Elevation (°)"] = np.round(altitudes, 2)
        merged["Slant Range (km)"] = np.round(slant_ranges, 2)
        return merged

    def _add_pfd(self, frame: pd.DataFrame) -> None:
        link = self.link
        if link.pfd_limit_low is None or not link.occupied_bandwidth_hz or "EIRP (dBW)" not in frame:
            return
        slant_m = frame["Slant Range (km)"].to_numpy(dtype=float) * 1000.0
        eirp = frame["EIRP (dBW)"].to_numpy(dtype=float)
        pfd = eirp - 10 * np.log10(4 * np.pi) - 20 * np.log10(slant_m) - 10 * np.log10(
            link.occupied_bandwidth_hz / 4000.0
        )
        limit = calculations.pfd_limit_at_elevation(
            link.pfd_limit_low, link.pfd_limit_high, frame["Elevation (°)"].to_numpy(dtype=float)
        )
        frame["PFD (dBW/m²/4kHz)"] = np.round(pfd, 2)
        frame["PFD Limit (dBW/m²/4kHz)"] = np.round(limit, 2)
        frame["PFD Margin (dB)"] = np.round(limit - pfd, 2)

    def contacts(self, frame: pd.DataFrame) -> list[tuple[datetime, datetime, str]]:
        """Contact windows (start, end, A/D) of a day of samples.

        The end is the first sample out of contact (or the last sample of the
        day for a contact still open at midnight).
        """
        times = list(frame["Time (UTC)"])
        visible = frame["Visible"].to_numpy(dtype=bool)
        windows = []
        start = None
        for i, now in enumerate(visible):
            if now and start is None:
                start = times[i]
            elif not now and start is not None:
                windows.append((start, times[i], pass_direction(self.sat, self.ts, start, times[i])))
                start = None
        if start is not None:
            windows.append((start, times[-1], pass_direction(self.sat, self.ts, start, times[-1])))
        return windows


def utc_midnight(date) -> datetime:
    return datetime(date.year, date.month, date.day, tzinfo=timezone.utc)


def _number(payload: dict, key: str, default=None):
    value = payload.get(key, default)
    return default if value in (None, "") else float(value)


def link_inputs_from_parameters(params: dict) -> LinkInputs:
    """The link inputs of a mission's parameters.json (as the window loads them).

    The static link budget sections give the losses, the Tx chain, the PFD
    limits and the required Eb/No, as with the "Link Losses from the Static
    Link Budget" options of the Downlink and Uplink tabs switched on.
    """
    freq = _number(params, "frequency_ghz")
    d_gs = _number(params, "antenna_diameter_m", 3.0)
    dl_section = params.get(calculations.FIXED_LINK_BUDGET_SECTIONS["downlink"]) or {}
    ul_section = params.get(calculations.FIXED_LINK_BUDGET_SECTIONS["uplink"]) or {}
    dl_values = static_link_values(dl_section, freq, d_gs) if dl_section else None
    uplink = None
    ul_freq = _number(params, "uplink_frequency_ghz")
    if ul_freq is not None and params.get("eirp_gs_dbw") is not None and params.get("gt_sat_dbk") is not None:
        ul_bitrate = _number(params, "uplink_bitrate_mbps")
        uplink = UplinkInputs(
            freq_ghz=ul_freq,
            eirp_gs=_number(params, "eirp_gs_dbw"),
            gt_sat=_number(params, "gt_sat_dbk"),
            demod_loss=_number(params, "demod_loss_sat_db"),
            other_att=_number(params, "uplink_other_attenuations_db"),
            bitrate_bps=ul_bitrate * 1e6 if ul_bitrate is not None else None,
            overhead=_number(params, "uplink_overhead"),
        )
    ul_values = static_link_values(ul_section, ul_freq, d_gs) if (ul_section and ul_freq) else None
    tsys = _number(params, "tsys_ref_k")
    gt_model = None
    if tsys is not None:
        gt_model = {
            "ref_elevation_deg": _number(params, "gt_ref_elevation_deg", calculations.MIN_ELEVATION_DEG),
            "tsys_ref_k": tsys,
            "radome_loss_db": _number(params, "radome_loss_db", 0.0),
        }
    bitrate = _number(params, "bitrate_mbps") * 1e6
    rolloff = _number(params, "rolloff")
    spectral_eff = _number(params, "spectral_efficiency_bpshz")
    return LinkInputs(
        freq_ghz=freq,
        link_availability_pct=_number(params, "link_availability_pct"),
        eirp_sat=_number(params, "eirp_sat_dbw"),
        gt_gs=_number(params, "gt_gs_dbk"),
        demod_loss=_number(params, "demod_loss_db"),
        bitrate_bps=bitrate,
        overhead=_number(params, "overhead"),
        other_att=_number(params, "other_attenuations_db", 0.0),
        cio_db=_number(params, "c_io_dbhz"),
        scintillation=False,
        gt_ref_diameter_m=_number(params, "gt_ref_diameter_m"),
        gt_model=gt_model,
        loss_kwargs=pass_loss_kwargs(dl_values, ul_values if uplink else None),
        uplink=uplink,
        pfd_limit_low=dl_values["pfd_limit"] if dl_values else None,
        pfd_limit_high=dl_values["pfd_limit_high"] if dl_values else None,
        occupied_bandwidth_hz=bitrate * (1 + rolloff) / spectral_eff if rolloff is not None and spectral_eff else None,
        sun_noise=bool(params.get("sun_noise", True)),
        solar_flux_sfu=_number(params, "solar_flux_sfu"),
    )


def station_from_catalogue(name: str, default_diameter_m: float = 3.0) -> Station:
    """A station of ground_stations.txt, with its diameter and horizon mask file if any."""
    lat, lon, alt_m = calculations.GROUND_STATIONS[name]
    mask_path = calculations.horizon_mask_file(name)
    return Station(
        name, lat, lon, alt_m,
        d_gs=calculations.GROUND_STATION_DIAMETERS.get(name, default_diameter_m),
        mask=calculations.load_horizon_mask(mask_path) if mask_path else None,
    )


def constellation_tles(orbit: dict, epoch: datetime, name: str) -> list[tuple[str, tuple[str, str]]]:
    """Synthetic TLEs of every satellite of a custom-orbit constellation.

    ``orbit`` is the orbit section of the mission's parameters: the planes
    are its LTAN list; ``satellites_per_plane`` satellites (default 1) are
    spaced evenly in each plane, and ``phasing_deg`` (default 0) shifts the
    argument of latitude from one plane to the next. Returns (name, TLE
    lines) per satellite, e.g. "EPS-Sterna 15:30 #1".
    """
    import tle_source

    planes = orbit.get("ltan") or ["22:20"]
    planes = planes if isinstance(planes, list) else [planes]
    per_plane = max(1, int(orbit.get("satellites_per_plane") or 1))
    phasing = float(orbit.get("phasing_deg") or 0.0)
    result = []
    for p, plane in enumerate(planes):
        for j in range(per_plane):
            tle = tle_source.synthetic_tle(
                float(orbit.get("altitude_km", 628.0)),
                tle_source.ltan_hours(plane),
                epoch,
                orbit.get("inclination_deg"),
                name=f"{name} {plane} #{j + 1}".upper(),
                norad_id=99900 + p * 10 + j,
                eccentricity=float(orbit.get("eccentricity") or 0.0),
                arg_perigee_deg=float(orbit.get("arg_perigee_deg") or 0.0),
                arg_latitude_deg=(j * 360.0 / per_plane + p * phasing) % 360.0,
            )
            result.append((f"{name} {plane} #{j + 1}", (tle.line1, tle.line2)))
    return result
