# Constants and core calculations for the satellite link budget tool
import numpy as np
import itur as itu
import astropy.units as u
from scipy.interpolate import interp1d
from typing import Any, Dict, Optional

# Define ground stations as a constant dictionary
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


def antenna_pattern(angle_deg: float) -> float:
    """Return the antenna pointing loss for a given off-boresight angle."""
    gain_at_angle = PATTERN_INTERP(angle_deg)
    loss = MAX_ANTENNA_GAIN - gain_at_angle
    return max(loss, 0.0)


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
        Rain rate exceeded for 0.01\% of the time in mm/h.
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

    Returns
    -------
    dict
        Dictionary containing computed parameters with the following keys:

        ``"Time (UTC)"`` : :class:`datetime.datetime`
            Timestamp corresponding to ``t_sky``.
        ``"Elevation (°)"`` : float
            Elevation angle in degrees.
        ``"Slant Range (km)"`` : float or ``None``
            Distance to the satellite when visible.
        ``"Path Loss (dB)"`` : float or ``None``
            Free-space path loss.
        ``"Pointing Loss (dB)"`` : float or ``None``
            Loss due to antenna off-pointing.
        ``"Off Boresight Angle (°)"`` : float or ``None``
            Angle between satellite boresight and ground station direction.
        ``"Rx Power (dBW)"`` : float or ``None``
            Received carrier power.
        ``"C/(No+Io) (dBHz)"`` : float or ``None``
            Carrier-to-noise-plus-interference density.
        ``"Eb/No (dB)"`` : float or ``None``
            Energy-per-bit to noise density.
        ``"Visible"`` : str
            ``"YES"`` if the elevation is above 5°; otherwise ``"NO"``.
    """
    diff = sat - gs
    topocentric = diff.at(t_sky)
    alt, az, dist = topocentric.altaz()
    elev = alt.degrees
    visible = elev >= 5

    slant_range_km = path_loss = atm_att = None
    rx_power = cn0 = ebno = off_boresight_angle = pointing_loss = None

    if visible:
        slant_range_km = topocentric.distance().to(u.km).value
        path_loss = (
            20 * np.log10(slant_range_km)
            + 20 * np.log10(freq.to(u.GHz).value)
            + 92.45
        )

        try:
            Ag, Ac, Ar, As, A_tot = itu.atmospheric_attenuation_slant_path(
                gs.latitude.degrees,
                gs.longitude.degrees,
                freq.to(u.GHz).value,
                elev,
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
            atm_att = float(A_tot.value if hasattr(A_tot, 'value') else A_tot)
        except Exception:
            atm_att = 0.0

        r_sat = sat.at(t_sky).position.km
        r_gs = gs.at(t_sky).position.km
        bore = -r_sat / np.linalg.norm(r_sat)
        to_gs = (r_gs - r_sat) / np.linalg.norm(r_gs - r_sat)
        angle_rad = np.arccos(np.clip(np.dot(bore, to_gs), -1.0, 1.0))
        off_boresight_angle = np.degrees(angle_rad)
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
        "Pointing Loss (dB)": pointing_loss,
        "Off Boresight Angle (°)": off_boresight_angle,
        "Rx Power (dBW)": rx_power,
        "C/(No+Io) (dBHz)": cn0,
        "Eb/No (dB)": ebno,
        "Visible": "YES" if visible else "NO",
    }
