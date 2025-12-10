"""Constants and helper functions for the satellite link budget tool.

This module provides antenna pattern utilities, Doppler and atmospheric
attenuation calculations, and functions to compute detailed link budget
parameters.
"""
import numpy as np
import itur as itu
import sys
import astropy.units as u
import os
from scipy.interpolate import interp1d
from typing import Any, Dict, Optional

# Minimum elevation angle (degrees) for visibility calculations
MIN_ELEVATION_DEG = 5.0

# Define ground stations loaded from external file. When running from a PyInstaller
# bundle, prefer a file placed alongside the executable (or in the current working
# directory) so users can override the packaged defaults without rebuilding.
def _resolve_ground_stations_file() -> str:
    """Pick the most appropriate ground stations file path."""

    candidates = []

    # Highest priority: an explicit environment override.
    env_path = os.environ.get("GROUND_STATIONS_FILE")
    if env_path:
        candidates.append(env_path)

    # If bundled, allow a file next to the executable to override packaged data.
    if getattr(sys, "frozen", False):
        exe_dir = os.path.dirname(sys.executable)
        candidates.append(os.path.join(exe_dir, "ground_stations.txt"))

    # Local working directory comes next so users can drop in a custom file while
    # running from source.
    candidates.append(os.path.join(os.getcwd(), "ground_stations.txt"))

    # Fallback to the repository/module path.
    candidates.append(os.path.join(os.path.dirname(__file__), "ground_stations.txt"))

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate

    # If none exist, still return the last candidate so the caller can emit a
    # meaningful error.
    return candidates[-1]


# Define ground stations loaded from external file
GROUND_STATIONS_FILE = _resolve_ground_stations_file()


def load_ground_stations(file_path: str | None = None) -> dict[str, tuple[float, float, float]]:
    """Load ground stations from a text file.

    Each non-empty line must contain four comma-separated values:
    station name, latitude (deg), longitude (deg) and altitude (m).
    Lines starting with ``#`` are treated as comments and ignored.
    """

    # 🔴 NOVITÀ: se non viene passato nulla, usa il valore *corrente*
    # di GROUND_STATIONS_FILE (che la GUI può aggiornare a runtime).
    if file_path is None:
        file_path = GROUND_STATIONS_FILE

    stations: dict[str, tuple[float, float, float]] = {}
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Ground station file not found: {file_path}")

    with open(file_path, "r", encoding="utf-8") as f:
        for line_no, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) != 4:
                raise ValueError(
                    f"Line {line_no} in {file_path} must have 4 comma-separated fields"
                )
            name, lat_str, lon_str, alt_str = parts
            try:
                lat = float(lat_str)
                lon = float(lon_str)
                alt = float(alt_str)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid numeric value on line {line_no} in {file_path}: {exc}"
                ) from exc
            stations[name] = (lat, lon, alt)

    if not stations:
        raise ValueError(f"No ground stations found in file: {file_path}")
    return stations


try:
    GROUND_STATIONS = load_ground_stations()
except Exception as exc:  # pragma: no cover - fallback for packaging/runtime errors
    print(f"Warning: using fallback ground stations: {exc}", file=sys.stderr)
    GROUND_STATIONS = {
        "Darmstadt": (49.8700, 8.6500, 144),
        "Lannion": (48.7333, -3.4542, 31),
        "Maspalomas": (27.7614, -15.5865, 250),
        "Athens": (37.9838, 23.7275, 70),
        "Kangerlussuaq": (67.0121, -50.7078, 50),
        "Svalbard": (78.2232, 15.6267, 10),
    }


# Antenna pattern data as constants
ANTENNA_PATTERN_ANGLES = np.array([
    -65, -55, -45, -35, -25, -15, -5, 0, 5, 15, 25, 35, 45, 55, 65
])
ANTENNA_PATTERN_GAINS = np.array([
    0.5, 1.2, 2.3, 3.5, 4.3, 5.0, 5.5, 5.6, 5.5, 5.0, 4.3, 3.5, 2.3, 1.2, 0.5
])
PATTERN_INTERP = interp1d(ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS,
                          kind='linear', fill_value=0.0, bounds_error=False)
MAX_ANTENNA_GAIN = 5.6


def load_antenna_pattern(file_path: str) -> None:
    """Load antenna pattern data from a CSV or text file.

    The file must contain two numeric columns: angle in degrees and gain in dB.
    Existing global pattern arrays are replaced and the interpolation function
    is updated accordingly.
    """

    global ANTENNA_PATTERN_ANGLES, ANTENNA_PATTERN_GAINS, PATTERN_INTERP, MAX_ANTENNA_GAIN

    data = np.loadtxt(file_path, delimiter=",")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("Antenna pattern file must have at least two columns")

    ANTENNA_PATTERN_ANGLES = data[:, 0]
    ANTENNA_PATTERN_GAINS = data[:, 1]
    PATTERN_INTERP = interp1d(
        ANTENNA_PATTERN_ANGLES,
        ANTENNA_PATTERN_GAINS,
        kind="linear",
        fill_value=0.0,
        bounds_error=False,
    )
    MAX_ANTENNA_GAIN = float(np.max(ANTENNA_PATTERN_GAINS))
# Speed of light in km/s used for Doppler calculations
C_KM_S = 299_792.458
def antenna_pattern(angle_deg: float | np.ndarray) -> float | np.ndarray:
    """Return the antenna pointing loss for a given off-boresight angle.

    Parameters
    ----------
    angle_deg : float or ndarray
        Off-boresight angle(s) in degrees.

    Returns
    -------
    float or ndarray
        Corresponding pointing loss in dB. If an array of angles is given,
        an array of losses with the same shape is returned.
    """

    angles = np.atleast_1d(angle_deg)
    gains = PATTERN_INTERP(angles)
    losses = MAX_ANTENNA_GAIN - gains
    losses = np.clip(losses, 0.0, None)

    if np.isscalar(angle_deg):
        return float(losses)
    return losses



def calculate_doppler_shift(
    topocentric: "ToposAt",  # type: ignore
    freq: u.Quantity,
) -> float:
    """Compute Doppler shift in Hz for the given topocentric position.

    Parameters
    ----------
    topocentric : :class:`~skyfield.positionlib.Geocentric`
        Relative position of satellite with respect to ground station.
    freq : :class:`~astropy.units.Quantity`
        Transmit frequency.

    Returns
    -------
    float
        Doppler frequency shift in Hz. Positive values indicate an
        approaching satellite (frequency increase).
    """
    los = topocentric.position.km
    rel_vel = topocentric.velocity.km_per_s
    los_unit = los / np.linalg.norm(los)
    radial_velocity = np.dot(rel_vel, los_unit)
    doppler_hz = -radial_velocity / C_KM_S * freq.to(u.Hz).value
    doppler_khz = doppler_hz/1000
    return float(doppler_khz)

def atmospheric_attenuation(
    lat: float,
    lon: float,
    freq_ghz: float,
    p: float,
    d_gs: float,
    alt_gs: float,
    r001: float,
) -> float:
    """Return slant path atmospheric attenuation in dB."""
    try:
        _, _, _, _, A_tot = itu.atmospheric_attenuation_slant_path(
            lat,
            lon,
            freq_ghz,
            MIN_ELEVATION_DEG,
            p,
            d_gs,
            hs=alt_gs,
            R001=r001,
            return_contributions=True,
            include_gas=True,
            include_rain=True,
            include_clouds=True,
            include_scintillation=True,
        )
        return float(A_tot.value if hasattr(A_tot, "value") else A_tot)
    except Exception as e:
        print(
            f"Error calculating atmospheric attenuation: {e}",
            file=sys.stderr,
        )
        print(
            "  Inputs for error: "
            f"lat={lat}, lon={lon}, freq_GHz={freq_ghz}, "
            f"elev={MIN_ELEVATION_DEG}, P={p}, "
            f"R001={r001}, D={d_gs}, h_s={alt_gs}",
            file=sys.stderr,
        )
        return 0.0


def prepare_topocentric_data(
    sat: "Satellite",  # type: ignore
    gs: "GroundStation",  # type: ignore
    times: "Time",  # type: ignore
    freq: u.Quantity,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pre-compute geometry values for a sequence of times.

    Parameters
    ----------
    sat, gs : Skyfield objects
        Satellite and ground station used for the calculations.
    times : :class:`~skyfield.timelib.Time`
        Array of times for which to compute the parameters.
    freq : :class:`~astropy.units.Quantity`
        Downlink frequency used for Doppler computation.

    Returns
    -------
    tuple of ndarrays
        ``(altitudes_deg, slant_ranges_km, doppler_khz, off_boresight_deg)``.
    """

    diff = (sat - gs).at(times)
    alt, _, dist = diff.altaz()
    altitudes_deg = alt.degrees
    slant_ranges_km = dist.km

    los = diff.position.km
    rel_vel = diff.velocity.km_per_s
    los_unit = los / np.linalg.norm(los, axis=0)
    radial_velocity = np.sum(rel_vel * los_unit, axis=0)
    doppler_hz = -radial_velocity / C_KM_S * freq.to(u.Hz).value
    doppler_khz = doppler_hz / 1000.0

    r_sat = sat.at(times).position.km
    r_gs = gs.at(times).position.km
    bore = -r_sat / np.linalg.norm(r_sat, axis=0)
    to_gs = r_gs - r_sat
    to_gs_unit = to_gs / np.linalg.norm(to_gs, axis=0)
    cos_angle = np.sum(bore * to_gs_unit, axis=0)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    off_boresight_deg = np.degrees(np.arccos(cos_angle))

    return altitudes_deg, slant_ranges_km, doppler_khz, off_boresight_deg




def calculate_link_budget_parameters(
t_sky: "Time", # type: ignore
    sat: "Satellite", # type: ignore
    gs: "GroundStation", # type: ignore
    freq: u.Quantity,
    p: float,
    r001: float,
    d_gs: float,
    alt_gs: float,
    eirp_sat: float,
    gt_gs: float,
    demod_loss: float,
    bitrate: float,
    overhead: float,
    cisat_lin: Optional[float],
    other_att: float,
    atm_att: Optional[float] = None,
    pre_alt_deg: Optional[float] = None,
    pre_slant_range_km: Optional[float] = None,
    pre_doppler_khz: Optional[float] = None,
    pre_off_boresight_angle: Optional[float] = None,
) -> Dict[str, Any]:
    
    """Compute link budget parameters for a single time step.
    Parameters
    ----------
    t_sky : :class:`~skyfield.timelib.Time`
        Skyfield time object for which the link budget is evaluated.
    sat : :class:`~skyfield.api.EarthSatellite`
        Satellite for which the pass is computed.
    gs : :class:`~skyfield.api.wgs84.GeographicPosition`
        Ground station location.
    freq : :class:`~astropy.units.Quantity`
        Downlink frequency.
    p : float
        Percentage of time attenuation is exceeded (100 - link availability).
    r001 : float
        Rain rate exceeded for 0.01% of the time in mm/h.
    d_gs : float
        Ground station antenna diameter in metres.
    alt_gs : float
        Ground station altitude in kilometres.
    eirp_sat : float
        Satellite EIRP in dBW.
    gt_gs : float
        Ground station G/T in dB/K.
    demod_loss : float
        Demodulator implementation loss in dB.
    bitrate : float
        Channel bit rate in bits per second.
    overhead : float
        Coding overhead factor used for Eb/No calculation.
    cisat_lin : float or None
        Satellite C/I as a linear ratio, or ``None`` if not used.
    other_att : float
        Any other static attenuation to subtract from the link budget in dB.
    atm_att : float or None
        Pre-computed atmospheric attenuation in dB. If ``None``, it will be
        calculated for each call.
    pre_alt_deg, pre_slant_range_km, pre_doppler_khz, pre_off_boresight_angle : float or None
        Optional pre-computed geometry parameters for this time step. Providing
        these values avoids repeated calls to Skyfield when iterating over many
        epochs.

    Returns
    -------
    dict
        Dictionary containing computed parameters with the following keys:

        ``"Time (UTC)"`` : :class:`datetime.datetime`
            Timestamp corresponding to ``t_sky``.
        ``"Elevation (°)"`` : float
            Elevation angle in degrees.
        ``"Slant Range (km)"`` : float
            Distance to the satellite.
        ``"Path Loss (dB)"`` : float
            Free-space path loss.
         ``"Pointing Loss (dB)"`` : float
            Loss due to antenna off-pointing.
        ``"Off Boresight Angle (°)"``
            Angle between satellite boresight and ground station direction.
        ``"Rx Power (dBW)"`` : float
            Received carrier power.
        ``"C/(No+Io) (dBHz)"`` : float
            Carrier-to-noise-plus-interference density.
        ``"Eb/No (dB)"`` : float
            Energy-per-bit to noise density.
        ``"Doppler Shift (Hz)"`` : float
            Instantaneous Doppler frequency shift.            
        ``"Visible"`` : str
            ``"YES"`` if the elevation is above ``MIN_ELEVATION_DEG``°; otherwise ``"NO"``.
    """
    if pre_alt_deg is None or pre_slant_range_km is None or pre_doppler_khz is None or pre_off_boresight_angle is None:
        diff = sat - gs
        topocentric = diff.at(t_sky)
    if pre_doppler_khz is None:
        doppler_khz = calculate_doppler_shift(topocentric, freq)
    else:
        doppler_khz = pre_doppler_khz
    if pre_alt_deg is None or pre_slant_range_km is None:
        alt, _, dist = topocentric.altaz()
        elev = alt.degrees
        slant_range_km = dist.km
    else:
        elev = pre_alt_deg
        slant_range_km = pre_slant_range_km
    visible = elev >= MIN_ELEVATION_DEG

    path_loss = (
        20 * np.log10(slant_range_km)
        + 20 * np.log10(freq.to(u.GHz).value)
        + 92.45
    )

    if atm_att is None:
        atm_att = atmospheric_attenuation(
            gs.latitude.degrees,
            gs.longitude.degrees,
            freq.to(u.GHz).value,
            p,
            d_gs,
            alt_gs,
            r001,
        )

    if pre_off_boresight_angle is None:
        r_sat = sat.at(t_sky).position.km
        r_gs = gs.at(t_sky).position.km
        bore = -r_sat / np.linalg.norm(r_sat)
        to_gs = (r_gs - r_sat) / np.linalg.norm(r_gs - r_sat)
        angle_rad = np.arccos(np.clip(np.dot(bore, to_gs), -1.0, 1.0))
        off_boresight_angle = np.degrees(angle_rad)
    else:
        off_boresight_angle = pre_off_boresight_angle
    pointing_loss = antenna_pattern(off_boresight_angle)

    rx_power = eirp_sat - path_loss - atm_att - other_att - pointing_loss
    cno_db = rx_power + gt_gs + 228.6 - demod_loss

    if cisat_lin is not None:
        cno_lin = 10 ** (cno_db / 10.0)
        cni_lin = 1.0 / (1.0 / cno_lin + 1.0 / cisat_lin)
        cn0 = 10.0 * np.log10(cni_lin)
    else:
        cn0 = cno_db

    if cn0 is not None and bitrate > 0 and overhead > 0:
        ebno = cn0 - 10 * np.log10(bitrate / overhead)
    else:
        ebno = np.nan

    return {
        "Time (UTC)": t_sky.utc_datetime(),
        "Elevation (°)": elev,
        "Slant Range (km)": slant_range_km,
        "Path Loss (dB)": path_loss,
        "Atmospheric Att (dB)": atm_att,
        "Pointing Loss (dB)": pointing_loss,
        "Off Boresight Angle (°)": off_boresight_angle,
        "Rx Power (dBW)": rx_power,
        "C/(No+Io) (dBHz)": cn0,
        "Eb/No (dB)": ebno,
        "Doppler Shift (kHz)": doppler_khz,        
        "Visible": "YES" if visible else "NO",
    }
