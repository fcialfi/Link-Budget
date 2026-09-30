"""Constants and helper functions for the satellite link budget tool.

This module provides antenna pattern utilities, Doppler and atmospheric
attenuation calculations, and functions to compute detailed link budget
parameters.
"""
import numpy as np
import itur as itu
import sys
import astropy.units as u
import json
import os
from scipy.interpolate import interp1d
from typing import Any, Dict, Optional

# Minimum elevation angle (degrees) for visibility calculations
MIN_ELEVATION_DEG = 5.0

# Mean equatorial Earth radius (km, WGS84), used for the TLE-independent
# geometric slant range calculation.
EARTH_RADIUS_KM = 6378.137

# Define ground stations loaded from an external text file. The file is **not**
# packaged into the PyInstaller bundle so it can be edited or replaced after
# compilation. When running a frozen executable, the loader first looks for a
# `ground_stations.txt` next to the EXE (or in the current working directory)
# and only falls back to the repository copy when running from source.
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

    # Fallback to the repository/module path only when running from source.
    if not getattr(sys, "frozen", False):
        candidates.append(os.path.join(os.path.dirname(__file__), "ground_stations.txt"))

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate

    # If none exist, still return the last candidate so the caller can emit a
    # meaningful error.
    return candidates[-1]


def resolve_optional_data_file(filename: str, env_var: str | None = None) -> str | None:
    """Find a companion data file, or return ``None`` if it isn't present.

    Looks in the same locations as the ground-station catalogue (an
    optional environment variable override, next to a frozen executable,
    the current working directory, and finally next to this module when
    running from source), but -- unlike the ground stations file -- there
    is no bundled fallback, so a missing file is a normal, silent case.
    """

    candidates = []
    if env_var:
        env_path = os.environ.get(env_var)
        if env_path:
            candidates.append(env_path)
    if getattr(sys, "frozen", False):
        exe_dir = os.path.dirname(sys.executable)
        candidates.append(os.path.join(exe_dir, filename))
    candidates.append(os.path.join(os.getcwd(), filename))
    if not getattr(sys, "frozen", False):
        candidates.append(os.path.join(os.path.dirname(__file__), filename))

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate
    return None


# Define ground stations loaded from external file
GROUND_STATIONS_FILE = _resolve_ground_stations_file()


def load_ground_stations(file_path: str | None = None) -> dict[str, tuple[float, float, float]]:
    """Load ground stations from a text file.

    Each non-empty line must contain four comma-separated values:
    station name, latitude (deg), longitude (deg) and altitude (m).
    Lines starting with ``#`` are treated as comments and ignored.
    """

    # Default to the current GROUND_STATIONS_FILE, which the GUI may have
    # updated at runtime via reload_ground_stations().
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


def reload_ground_stations(file_path: str) -> dict[str, tuple[float, float, float]]:
    """Reload the global ground station list from ``file_path``.

    This helper validates the file by reusing :func:`load_ground_stations`, then
    updates :data:`GROUND_STATIONS_FILE` and :data:`GROUND_STATIONS` in-place so
    other modules (e.g. the GUI) immediately see the refreshed catalogue.
    """

    global GROUND_STATIONS_FILE, GROUND_STATIONS

    file_path = os.path.abspath(file_path)
    stations = load_ground_stations(file_path)
    GROUND_STATIONS_FILE = file_path
    GROUND_STATIONS = stations
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
        return float(losses[0])
    return losses



def calculate_doppler_shift(
    topocentric: "ToposAt",  # type: ignore
    freq: u.Quantity,
) -> float:
    """Compute Doppler shift in kHz for the given topocentric position.

    Parameters
    ----------
    topocentric : :class:`~skyfield.positionlib.Geocentric`
        Relative position of satellite with respect to ground station.
    freq : :class:`~astropy.units.Quantity`
        Transmit frequency.

    Returns
    -------
    float
        Doppler frequency shift in **kHz**. Positive values indicate an
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
    include_scintillation: bool = True,
) -> float:
    """Return slant path atmospheric attenuation in dB.

    The ITU-R model internally estimates the point rainfall rate (R001) from
    Recommendation P.837 when it is not provided explicitly.

    ``include_scintillation`` toggles the ITU-R P.618 *tropospheric*
    scintillation term (relevant mainly above ~4 GHz). This is unrelated to
    ionospheric scintillation, which ``itur`` does not model.
    """
    try:
        _, _, _, _, A_tot = itu.atmospheric_attenuation_slant_path(
            lat,
            lon,
            freq_ghz,
            MIN_ELEVATION_DEG,
            p,
            d_gs,
            hs=alt_gs,
            return_contributions=True,
            include_gas=True,
            include_rain=True,
            include_clouds=True,
            include_scintillation=include_scintillation,
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
            f"D={d_gs}, h_s={alt_gs}",
            file=sys.stderr,
        )
        return 0.0


# Sections of ``parameters.json`` holding the inputs of the fixed-elevation
# (static) link budget popups, so they don't have to be retyped every time.
FIXED_LINK_BUDGET_SECTIONS = {
    "downlink": "fixed_elevation_downlink",
    "uplink": "fixed_elevation_uplink",
}


def extract_json_section(payload: Any, section: str) -> Dict[str, Any]:
    """Return ``payload[section]`` from a parameters JSON payload.

    Raises ``ValueError`` when the payload is not a JSON object or has no
    such section (or the section is not an object).
    """

    if not isinstance(payload, dict):
        raise ValueError("The file must contain a JSON object.")
    if section not in payload:
        raise ValueError(f"The file has no '{section}' section.")
    values = payload[section]
    if not isinstance(values, dict):
        raise ValueError(f"The '{section}' section must be a JSON object.")
    return values


def save_json_section(file_path: str, section: str, values: Dict[str, Any]) -> None:
    """Write ``values`` as ``section`` of the JSON file at ``file_path``.

    Every other key already in the file is kept, so the static link budget
    inputs can be saved straight into ``parameters.json`` alongside the main
    parameters (and are then auto-loaded at startup with it). A missing or
    empty file is created from scratch; an existing file that is not a JSON
    object raises ``ValueError`` instead of being overwritten.
    """

    payload: Dict[str, Any] = {}
    if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
        with open(file_path, "r", encoding="utf-8") as f:
            existing = json.load(f)
        if not isinstance(existing, dict):
            raise ValueError(f"'{file_path}' does not contain a JSON object; not overwriting it.")
        payload = existing
    payload[section] = values
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
        f.write("\n")


def slant_range_from_elevation(
    elevation_deg: float,
    sat_altitude_km: float,
    gs_altitude_km: float = 0.0,
) -> float:
    """Return the geometric slant range (km) for a fixed elevation angle.

    Assumes a spherical Earth and a circular orbit at ``sat_altitude_km``;
    does not require a TLE. Used by the fixed-elevation "preliminary" link
    budget calculators, where the actual pass geometry is deliberately
    abstracted away in favour of a single (typically worst-case, minimum)
    elevation angle -- the standard approach for a first-pass link margin
    check before running a full TLE-based analysis.
    """

    r_gs = EARTH_RADIUS_KM + gs_altitude_km
    r_sat = EARTH_RADIUS_KM + sat_altitude_km
    el_rad = np.radians(elevation_deg)
    return float(np.sqrt(r_sat**2 - (r_gs * np.cos(el_rad)) ** 2) - r_gs * np.sin(el_rad))


def atmospheric_attenuation_contributions(
    lat: float,
    lon: float,
    freq_ghz: float,
    elevation_deg: float,
    p: float,
    d_gs: float,
    alt_gs: float,
    include_rain: bool = True,
    include_clouds: bool = True,
    include_scintillation: bool = True,
) -> Dict[str, float]:
    """Return the individual ITU-R P.618 attenuation contributions, in dB.

    Unlike :func:`atmospheric_attenuation` (which always evaluates the model
    at ``MIN_ELEVATION_DEG`` as a single conservative margin applied across a
    whole TLE-derived pass), this accepts an explicit ``elevation_deg`` and
    breaks the total down into its gas/cloud/rain/scintillation components,
    for the fixed-elevation calculators that are independent of any TLE.
    """
    try:
        A_g, A_c, A_r, A_s, A_t = itu.atmospheric_attenuation_slant_path(
            lat,
            lon,
            freq_ghz,
            elevation_deg,
            p,
            d_gs,
            hs=alt_gs,
            return_contributions=True,
            include_gas=True,
            include_rain=include_rain,
            include_clouds=include_clouds,
            include_scintillation=include_scintillation,
        )

        def _val(x):
            return float(x.value if hasattr(x, "value") else x)

        return {
            "gas": _val(A_g),
            "cloud": _val(A_c),
            "rain": _val(A_r),
            "scintillation": _val(A_s),
            "total": _val(A_t),
        }
    except Exception as e:
        print(
            f"Error calculating atmospheric attenuation contributions: {e}",
            file=sys.stderr,
        )
        return {"gas": 0.0, "cloud": 0.0, "rain": 0.0, "scintillation": 0.0, "total": 0.0}


def vswr_mismatch_loss_db(vswr: float) -> float:
    """Return the transmitter/antenna mismatch (VSWR) loss in dB.

    Standard reflection-coefficient formula: with reflection coefficient
    ``gamma = (VSWR - 1) / (VSWR + 1)``, the mismatch loss is
    ``-10*log10(1 - gamma**2)``.
    """

    reflection = (vswr - 1.0) / (vswr + 1.0)
    return float(-10 * np.log10(1 - reflection**2))


def antenna_pointing_loss_db(depointing_deg: float, beamwidth_3db_deg: float) -> float:
    """Return the antenna pointing (depointing) loss in dB.

    Standard parabolic main-lobe approximation ``12 * (theta / theta_3dB)**2``,
    valid for depointing angles well inside the main lobe.
    """

    if beamwidth_3db_deg <= 0:
        raise ValueError("The 3 dB beamwidth must be positive.")
    return float(12.0 * (depointing_deg / beamwidth_3db_deg) ** 2)


def antenna_beamwidth_3db_deg(freq_ghz: float, diameter_m: float) -> float:
    """Approximate 3 dB beamwidth (deg) of a parabolic dish: ``70 * lambda / D``."""

    wavelength_m = 0.299792458 / freq_ghz
    return float(70.0 * wavelength_m / diameter_m)


def polarisation_mismatch_loss_db(tx_axial_ratio_db: float, rx_axial_ratio_db: float) -> Dict[str, float]:
    """Polarisation mismatch loss between two co-rotating elliptical antennas.

    Uses the classic polarisation-efficiency formula for two elliptically
    polarised antennas of the same sense, with the voltage axial ratios
    ``r = 10**(AR_dB / 20)``. The loss depends on the (unknown) relative
    orientation of the polarisation ellipses, so both the ``best`` (aligned
    major axes) and ``worst`` (orthogonal major axes) cases are returned.
    """

    r1 = 10 ** (tx_axial_ratio_db / 20.0)
    r2 = 10 ** (rx_axial_ratio_db / 20.0)
    denom = 2 * (1 + r1**2) * (1 + r2**2)
    cross = (r1**2 - 1) * (r2**2 - 1)
    best = 0.5 + (4 * r1 * r2 + cross) / denom
    worst = 0.5 + (4 * r1 * r2 - cross) / denom
    return {"best": float(-10 * np.log10(best)), "worst": float(-10 * np.log10(worst))}


TOLERANCE_DISTRIBUTIONS = ("TRI", "UNI", "GAU")


def tolerance_statistics(favourable_db: float, adverse_db: float, distribution: str) -> Dict[str, float]:
    """Mean offset and variance of a toleranced link budget contribution.

    ``favourable_db`` and ``adverse_db`` are the (non-negative) deviations of
    the parameter, expressed as their effect on the link margin, from its
    nominal value: the margin ranges over ``[-adverse_db, +favourable_db]``
    around nominal. ``distribution`` follows the usual ECSS/CCSDS link budget
    conventions:

    * ``TRI`` -- triangular, peaking at the nominal value;
    * ``UNI`` -- uniform between the adverse and favourable limits;
    * ``GAU`` -- Gaussian, with the adverse/favourable limits taken as
      +/-3 sigma.

    Returns ``mean_db`` (the offset of the mean from nominal) and
    ``variance_db2``.
    """

    a = -abs(adverse_db)
    b = abs(favourable_db)
    dist = distribution.upper()
    if dist == "TRI":
        mean = (a + b) / 3.0
        variance = (a**2 + b**2 - a * b) / 18.0
    elif dist == "UNI":
        mean = (a + b) / 2.0
        variance = (b - a) ** 2 / 12.0
    elif dist == "GAU":
        mean = (a + b) / 2.0
        variance = ((b - a) / 6.0) ** 2
    else:
        raise ValueError(f"Unknown tolerance distribution '{distribution}' (use TRI, UNI or GAU).")
    return {"mean_db": float(mean), "variance_db2": float(variance)}


def margin_statistics(nominal_margin_db: float, tolerances: list) -> Dict[str, Any]:
    """Statistical margin figures from a list of toleranced contributions.

    Each item of ``tolerances`` is a dict with ``name``, ``favourable_db``,
    ``adverse_db`` and ``distribution`` (see :func:`tolerance_statistics`).
    The contributions are assumed independent, so means add and variances
    add. Returns the favourable/adverse (all tolerances at their limits)
    margins, the mean and variance of the margin, ``mean_minus_3sigma_db``
    and ``worst_case_rss_db`` (nominal margin minus the root-sum-square of
    the adverse tolerances), plus the per-contribution breakdown under
    ``contributions``.
    """

    contributions = []
    mean_offset = 0.0
    variance = 0.0
    adverse_sq = 0.0
    fav_sum = 0.0
    adv_sum = 0.0
    for tol in tolerances:
        fav = abs(tol.get("favourable_db", 0.0))
        adv = abs(tol.get("adverse_db", 0.0))
        dist = tol.get("distribution", "TRI")
        stats = tolerance_statistics(fav, adv, dist)
        contributions.append(
            {
                "name": tol.get("name", ""),
                "favourable_db": fav,
                "adverse_db": adv,
                "distribution": dist.upper(),
                **stats,
            }
        )
        mean_offset += stats["mean_db"]
        variance += stats["variance_db2"]
        adverse_sq += adv**2
        fav_sum += fav
        adv_sum += adv

    mean = nominal_margin_db + mean_offset
    sigma = float(np.sqrt(variance))
    return {
        "nominal_db": float(nominal_margin_db),
        "favourable_db": float(nominal_margin_db + fav_sum),
        "adverse_db": float(nominal_margin_db - adv_sum),
        "mean_db": float(mean),
        "variance_db2": float(variance),
        "sigma_db": sigma,
        "mean_minus_3sigma_db": float(mean - 3 * sigma),
        "worst_case_rss_db": float(nominal_margin_db - np.sqrt(adverse_sq)),
        "contributions": contributions,
    }


def power_flux_density(
    eirp_dbw: float,
    slant_range_km: float,
    occupied_bandwidth_hz: Optional[float] = None,
) -> Dict[str, float]:
    """Power flux density (PFD) at the receiving site from EIRP and range alone.

    This is the geometric PFD (no atmospheric loss credit taken), matching
    the conservative convention used to check against regulatory PFD limits
    such as the ones in ECSS-E-ST-50-05C.

    Returns ``pfd_dbw_m2`` and, when ``occupied_bandwidth_hz`` is given,
    ``pfd_dbw_m2_per_4khz`` -- the PFD re-normalised to a 4 kHz reference
    bandwidth, assuming a uniform spectral density across the occupied
    bandwidth.
    """

    slant_range_m = slant_range_km * 1000.0
    pfd_dbw_m2 = eirp_dbw - 10 * np.log10(4 * np.pi) - 20 * np.log10(slant_range_m)
    result = {"pfd_dbw_m2": float(pfd_dbw_m2)}
    if occupied_bandwidth_hz:
        result["pfd_dbw_m2_per_4khz"] = float(pfd_dbw_m2 - 10 * np.log10(occupied_bandwidth_hz / 4000.0))
    return result


def calculate_fixed_elevation_link_budget(
    freq: u.Quantity,
    elevation_deg: float,
    sat_altitude_km: float,
    lat_gs: float,
    lon_gs: float,
    alt_gs_km: float,
    d_gs: float,
    eirp: Optional[float],
    gt: float,
    demod_loss: float,
    bitrate: float,
    overhead: float,
    other_att: float,
    pointing_loss_db: float,
    link_availability_pct: float,
    include_scintillation: bool = True,
    required_ebno: Optional[float] = None,
    tx_power_w: Optional[float] = None,
    antenna_circuit_loss_db: float = 0.0,
    vswr: Optional[float] = None,
    antenna_gain_dbi: Optional[float] = None,
    ionospheric_loss_db: float = 0.0,
    polarisation_loss_db: Optional[float] = 0.0,
    multipath_loss_db: float = 0.0,
    modulation_degradation_db: float = 0.0,
    occupied_bandwidth_hz: Optional[float] = None,
    pfd_limit_dbw_m2_4khz: Optional[float] = None,
    tx_axial_ratio_db: Optional[float] = None,
    rx_axial_ratio_db: Optional[float] = None,
    rx_depointing_deg: Optional[float] = None,
    rx_beamwidth_3db_deg: Optional[float] = None,
    formatting_overhead: float = 1.0,
    tolerances: Optional[list] = None,
) -> Dict[str, Any]:
    """Preliminary link budget at a single, user-chosen fixed elevation angle.

    This does not require a TLE: the slant range is derived purely from
    ``elevation_deg`` and ``sat_altitude_km`` via
    :func:`slant_range_from_elevation`. This is the standard "worst case"
    check used to verify a link margin exists before running a full
    TLE-based pass analysis, typically at the minimum operational elevation.

    The budget is computed twice: once for a clear-sky reference condition
    (gas absorption only -- rain, clouds and scintillation set to zero, the
    conventional "Clear" column in a preliminary link budget) and once for
    the rain-faded condition at the configured link availability
    (``p = 100 - link_availability_pct``, with rain/clouds/scintillation
    included). ``ionospheric_loss_db``, ``polarisation_loss_db`` and
    ``multipath_loss_db`` are static extra losses (independent of the
    clear/rain condition) subtracted from the received power in both
    columns, and ``modulation_degradation_db`` is added to the implementation
    loss alongside ``demod_loss``.

    When ``tx_power_w`` is given, the EIRP is computed from the transmit chain
    (``tx_power_dbw - antenna_circuit_loss_db - vswr_loss_db + antenna_gain_dbi``)
    and used in place of ``eirp`` for the whole budget (``eirp_source`` is
    then ``"tx_chain"``, otherwise ``"input"``); ``eirp`` may be ``None`` in
    that case. The EIRP actually used is returned as ``eirp_dbw``. When ``occupied_bandwidth_hz`` is given, the power
    flux density at the receiving site is also returned (see
    :func:`power_flux_density`), together with a margin against
    ``pfd_limit_dbw_m2_4khz`` when that is provided too.

    ``polarisation_loss_db=None`` derives the polarisation mismatch loss from
    ``tx_axial_ratio_db`` and ``rx_axial_ratio_db`` (worst-case orientation,
    see :func:`polarisation_mismatch_loss_db`) when both are given, and 0
    otherwise. ``rx_depointing_deg`` adds a receive antenna pointing loss
    (:func:`antenna_pointing_loss_db`), using ``rx_beamwidth_3db_deg`` or, if
    that is omitted, the beamwidth of a ``d_gs`` dish at ``freq``. Eb/No is
    referred to the bit rate including formatting but excluding coding,
    ``bitrate / overhead * formatting_overhead`` (``overhead`` being the
    total coded-over-information rate ratio). When ``tolerances`` (see
    :func:`margin_statistics`) and ``required_ebno`` are given, each
    condition also carries the statistical margin figures under
    ``margin_statistics``.

    Returns
    -------
    dict
        ``{"slant_range_km", "path_loss_db", "clear": {...}, "rain_faded": {...}, ...}``
        where each condition dict has ``gas_attenuation_db``,
        ``cloud_attenuation_db``, ``rain_attenuation_db``,
        ``scintillation_db``, ``atmospheric_attenuation_db``,
        ``total_propagation_loss_db``, ``rx_power_dbw``, ``cno_dbhz``,
        ``ebno_db`` and, when ``required_ebno`` is given, ``margin_db``. The
        top-level dict always carries ``ionospheric_loss_db``,
        ``polarisation_loss_db``, ``multipath_loss_db`` and
        ``modulation_degradation_db``, and conditionally
        ``tx_power_dbw``, ``vswr_loss_db``, ``effective_gain_dbi``,
        ``pfd_dbw_m2``, ``pfd_dbw_m2_per_4khz`` and ``pfd_margin_db``.
    """

    vswr_loss = vswr_mismatch_loss_db(vswr) if vswr else None
    effective_gain = None
    if antenna_gain_dbi is not None:
        effective_gain = float(antenna_gain_dbi - antenna_circuit_loss_db - (vswr_loss or 0.0))
    tx_power_dbw = None
    if tx_power_w is not None:
        if tx_power_w <= 0:
            raise ValueError("The transmitter power must be greater than 0 W.")
        if antenna_gain_dbi is None:
            raise ValueError("The antenna gain is required to compute the EIRP from the transmitter power.")
        tx_power_dbw = float(10 * np.log10(tx_power_w))
        eirp = tx_power_dbw + effective_gain
        eirp_source = "tx_chain"
    elif eirp is None:
        raise ValueError("Either the EIRP or the transmitter power (with antenna gain) is required.")
    else:
        eirp_source = "input"

    slant_range_km = slant_range_from_elevation(elevation_deg, sat_altitude_km, alt_gs_km)
    freq_ghz = freq.to(u.GHz).value
    path_loss = 20 * np.log10(slant_range_km) + 20 * np.log10(freq_ghz) + 92.45
    polarisation_range = None
    if tx_axial_ratio_db is not None and rx_axial_ratio_db is not None:
        polarisation_range = polarisation_mismatch_loss_db(tx_axial_ratio_db, rx_axial_ratio_db)
    if polarisation_loss_db is None:
        polarisation_loss_db = polarisation_range["worst"] if polarisation_range else 0.0

    rx_pointing_loss_db = 0.0
    if rx_depointing_deg is not None:
        if rx_beamwidth_3db_deg is None:
            rx_beamwidth_3db_deg = antenna_beamwidth_3db_deg(freq_ghz, d_gs)
        rx_pointing_loss_db = antenna_pointing_loss_db(rx_depointing_deg, rx_beamwidth_3db_deg)

    extra_static_loss_db = (
        ionospheric_loss_db + polarisation_loss_db + multipath_loss_db + rx_pointing_loss_db
    )
    demod_loss_total = demod_loss + modulation_degradation_db
    if bitrate > 0 and overhead > 0:
        ebno_bitrate = bitrate / overhead * (formatting_overhead or 1.0)
        bitrate_dbhz = float(10 * np.log10(ebno_bitrate))
    else:
        ebno_bitrate = 0.0
        bitrate_dbhz = float("nan")

    def _condition(p: float, include_rain: bool, include_clouds: bool, include_scint: bool) -> Dict[str, float]:
        atm = atmospheric_attenuation_contributions(
            lat_gs,
            lon_gs,
            freq_ghz,
            elevation_deg,
            p,
            d_gs,
            alt_gs_km,
            include_rain=include_rain,
            include_clouds=include_clouds,
            include_scintillation=include_scint,
        )
        rx_power = eirp - path_loss - atm["total"] - other_att - pointing_loss_db - extra_static_loss_db
        received_cno_db = rx_power + gt + 228.6
        cno_db = received_cno_db - demod_loss_total
        ebno = cno_db - bitrate_dbhz
        result: Dict[str, float] = {
            "gas_attenuation_db": atm["gas"],
            "cloud_attenuation_db": atm["cloud"],
            "rain_attenuation_db": atm["rain"],
            "scintillation_db": atm["scintillation"],
            "atmospheric_attenuation_db": atm["total"],
            "total_propagation_loss_db": path_loss + atm["total"] + ionospheric_loss_db + polarisation_loss_db,
            "rx_power_dbw": rx_power,
            "received_cno_dbhz": received_cno_db,
            "cno_dbhz": cno_db,
            "ebno_db": ebno,
        }
        if required_ebno is not None:
            result["margin_db"] = ebno - required_ebno
            if tolerances:
                result["margin_statistics"] = margin_statistics(result["margin_db"], tolerances)
        return result

    rain_p = max(0.001, min(50.0, 100.0 - link_availability_pct))
    out: Dict[str, Any] = {
        "slant_range_km": slant_range_km,
        "path_loss_db": path_loss,
        "clear": _condition(rain_p, include_rain=False, include_clouds=False, include_scint=False),
        "rain_faded": _condition(
            rain_p, include_rain=True, include_clouds=True, include_scint=include_scintillation
        ),
        "ionospheric_loss_db": ionospheric_loss_db,
        "polarisation_loss_db": polarisation_loss_db,
        "multipath_loss_db": multipath_loss_db,
        "modulation_degradation_db": modulation_degradation_db,
        "receiver_degradation_db": demod_loss,
        "rx_pointing_loss_db": rx_pointing_loss_db,
        "eirp_dbw": float(eirp),
        "eirp_source": eirp_source,
        "gt_dbk": gt,
        "frequency_mhz": freq_ghz * 1000.0,
        "bitrate_formatted_bps": ebno_bitrate,
        "bitrate_dbhz": bitrate_dbhz,
    }
    if polarisation_range is not None:
        out["polarisation_loss_range_db"] = polarisation_range
    if rx_depointing_deg is not None:
        out["rx_depointing_deg"] = rx_depointing_deg
        out["rx_beamwidth_3db_deg"] = rx_beamwidth_3db_deg
    if required_ebno is not None:
        out["required_ebno_db"] = required_ebno

    if vswr_loss is not None:
        out["vswr_loss_db"] = vswr_loss
    if effective_gain is not None:
        out["effective_gain_dbi"] = effective_gain
    if tx_power_dbw is not None:
        out["tx_power_dbw"] = tx_power_dbw

    pfd = power_flux_density(eirp, slant_range_km, occupied_bandwidth_hz)
    out["pfd_dbw_m2"] = pfd["pfd_dbw_m2"]
    out["pfd_dbm_m2"] = pfd["pfd_dbw_m2"] + 30.0
    if pfd_limit_dbw_m2_4khz is not None:
        out["pfd_limit_dbw_m2_4khz"] = pfd_limit_dbw_m2_4khz
    if "pfd_dbw_m2_per_4khz" in pfd:
        out["pfd_dbw_m2_per_4khz"] = pfd["pfd_dbw_m2_per_4khz"]
        if pfd_limit_dbw_m2_4khz is not None:
            out["pfd_margin_db"] = float(pfd_limit_dbw_m2_4khz - pfd["pfd_dbw_m2_per_4khz"])

    return out


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
    t_sky: "Time",  # type: ignore
    sat: "Satellite",  # type: ignore
    gs: "GroundStation",  # type: ignore
    freq: u.Quantity,
    p: float,
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
    uplink_freq: Optional[u.Quantity] = None,
    eirp_gs: Optional[float] = None,
    gt_sat: Optional[float] = None,
    demod_loss_ul: float | None = None,
    other_att_ul: float | None = None,
    atm_att_ul: Optional[float] = None,
    uplink_bitrate: float | None = None,
    uplink_overhead: float | None = None,
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
    uplink_freq : :class:`~astropy.units.Quantity`, optional
        Uplink frequency used to compute the reverse link budget. If omitted,
        uplink metrics are not calculated.
    eirp_gs : float, optional
        Ground station EIRP in dBW used for the uplink budget.
    gt_sat : float, optional
        Satellite G/T in dB/K used for the uplink budget.
    demod_loss_ul : float, optional
        Implementation loss on the satellite receiver for the uplink, in dB.
    other_att_ul : float, optional
        Additional static uplink attenuations in dB.
    atm_att_ul : float, optional
        Pre-computed uplink atmospheric attenuation in dB. If ``None`` it is
        recomputed using :func:`atmospheric_attenuation` with ``uplink_freq``.
    uplink_bitrate : float, optional
        Optional uplink bit rate in bits per second used to derive uplink
        Eb/No. If omitted, the downlink ``bitrate`` is reused.
    uplink_overhead : float, optional
        Optional uplink coding overhead. If omitted, the downlink ``overhead``
        factor is reused for uplink Eb/No calculations.

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
        ``"Doppler Shift (kHz)"`` : float
            Instantaneous Doppler frequency shift expressed in kilohertz.
        ``"Visible"`` : str
            ``"YES"`` if the elevation is above ``MIN_ELEVATION_DEG``°; otherwise ``"NO"``.

        When all uplink parameters are provided, the dictionary also includes
        ``"UL Path Loss (dB)"``, ``"UL Atmospheric Att (dB)"``,
        ``"UL Rx Power (dBW)"``, ``"UL C/No (dBHz)"`` and ``"UL Eb/No (dB)"``.
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

    results: Dict[str, Any] = {
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
        "UL Path Loss (dB)": None,
        "UL Atmospheric Att (dB)": None,
        "UL Pointing Loss (dB)": None,
        "UL Rx Power (dBW)": None,
        "UL C/No (dBHz)": None,
        "UL Eb/No (dB)": None,
    }

    if uplink_freq is not None and eirp_gs is not None and gt_sat is not None:
        ul_bitrate = bitrate if uplink_bitrate is None else uplink_bitrate
        ul_overhead = overhead if uplink_overhead is None else uplink_overhead

        ul_path_loss = (
            20 * np.log10(slant_range_km)
            + 20 * np.log10(uplink_freq.to(u.GHz).value)
            + 92.45
        )

        if atm_att_ul is None:
            atm_att_ul = atmospheric_attenuation(
                gs.latitude.degrees,
                gs.longitude.degrees,
                uplink_freq.to(u.GHz).value,
                p,
                d_gs,
                alt_gs,
            )

        ul_pointing_loss = pointing_loss
        ul_other_att = 0.0 if other_att_ul is None else other_att_ul
        ul_demod_loss = 0.0 if demod_loss_ul is None else demod_loss_ul

        ul_rx_power = eirp_gs - ul_path_loss - atm_att_ul - ul_other_att - ul_pointing_loss
        ul_cno = ul_rx_power + gt_sat + 228.6 - ul_demod_loss
        ul_ebno = ul_cno - 10 * np.log10(ul_bitrate / ul_overhead) if ul_bitrate > 0 and ul_overhead > 0 else np.nan

        results.update(
            {
                "UL Path Loss (dB)": ul_path_loss,
                "UL Atmospheric Att (dB)": atm_att_ul,
                "UL Pointing Loss (dB)": ul_pointing_loss,
                "UL Rx Power (dBW)": ul_rx_power,
                "UL C/No (dBHz)": ul_cno,
                "UL Eb/No (dB)": ul_ebno,
            }
        )

    return results
