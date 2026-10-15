"""Description of the mission parameters (``missions/<name>/parameters.json``).

Gives every known parameter a group, a label, a unit and a description, so
the window can show the file as a table and edit it without touching the
JSON. Parameters not in the schema are listed under "Other" with their key.
"""
from __future__ import annotations

from typing import Any, Optional

# (group, path in the JSON, label, unit, kind, description); kind: float, int, bool, str, list
SCHEMA: list[tuple[str, tuple[str, ...], str, str, str, str]] = [
    ("Mission", ("project_name",), "Project name", "", "str", "Shown in the window title and in the reports."),
    ("Mission", ("norad_id",), "NORAD ID", "", "int", "Catalogue number for the TLE download (launched satellites)."),
    ("Orbit", ("orbit", "altitude_km"), "Mean altitude", "km", "float", "Mean semi-major axis minus 6378.137 km."),
    ("Orbit", ("orbit", "eccentricity"), "Eccentricity", "", "float", "Mean eccentricity."),
    ("Orbit", ("orbit", "arg_perigee_deg"), "Argument of perigee", "deg", "float", "90 for a frozen orbit."),
    ("Orbit", ("orbit", "inclination_deg"), "Inclination", "deg", "float", "Empty: sun-synchronous."),
    ("Orbit", ("orbit", "ltan"), "Planes (LTAN)", "hh:mm", "list", "Local time of the ascending node of each plane."),
    ("Orbit", ("orbit", "satellites_per_plane"), "Satellites per plane", "", "int", "Evenly spaced in the plane."),
    ("Orbit", ("orbit", "phasing_deg"), "Phasing between planes", "deg", "float",
     "Shift of the argument of latitude from one plane to the next."),
    ("Orbit", ("orbit", "extra_satellites"), "Extra satellites", "LTAN@deg", "list",
     "Single satellites added to the regular ones, e.g. 19:30@120 (plane, argument of latitude at epoch)."),
    ("Downlink", ("frequency_ghz",), "Frequency", "GHz", "float", "Downlink carrier frequency."),
    ("Downlink", ("eirp_sat_dbw",), "EIRP SAT", "dBW", "float", "Satellite EIRP (when not from the Tx chain)."),
    ("Downlink", ("bitrate_mbps",), "Coded bit rate", "Mbps", "float", "Channel bit rate, coding included."),
    ("Downlink", ("overhead",), "Overhead (coded / info)", "", "float", "Coded bit rate / information bit rate."),
    ("Downlink", ("rolloff",), "Roll-off factor", "", "float", "Occupied bandwidth = coded rate × (1 + roll-off) / spectral efficiency."),
    ("Downlink", ("spectral_efficiency_bpshz",), "Spectral efficiency", "b/s/Hz", "float", "2 for QPSK."),
    ("Downlink", ("demod_loss_db",), "Demodulator loss", "dB", "float", "Receiver implementation loss."),
    ("Downlink", ("link_availability_pct",), "Link availability", "%", "float", "Rain and clouds for p = 100 − availability."),
    ("Downlink", ("other_attenuations_db",), "Other attenuations", "dB", "float", "Any further fixed loss."),
    ("Downlink", ("c_io_dbhz",), "C/Io", "dBHz", "float", "Carrier to interference density; empty: no interferer."),
    ("Ground station", ("gt_gs_dbk",), "G/T GS", "dB/K", "float", "G/T of the reference antenna, reference elevation, clear sky."),
    ("Ground station", ("gt_ref_diameter_m",), "G/T reference diameter", "m", "float", "Antenna the G/T GS refers to."),
    ("Ground station", ("gt_ref_elevation_deg",), "G/T reference elevation", "deg", "float", "Elevation of the G/T GS."),
    ("Ground station", ("tsys_ref_k",), "System noise temperature", "K", "float", "At the G/T reference elevation, clear sky."),
    ("Ground station", ("radome_loss_db",), "Radome loss", "dB", "float", "Subtracted from the G/T."),
    ("Ground station", ("sun_noise",), "Sun noise", "", "bool", "Sun in the antenna beam along the passes; empty: on."),
    ("Ground station", ("solar_flux_sfu",), "Solar flux", "sfu", "float", "For the Sun noise; empty: quiet Sun."),
    ("Ground station", ("antenna_diameter_m",), "Antenna diameter (default)", "m", "float",
     "For a station without a diameter in the station catalogue."),
    ("Ground station", ("min_elevation_deg",), "Minimum elevation", "deg", "float", "Minimum elevation of a contact."),
    ("Uplink", ("uplink_frequency_ghz",), "Frequency", "GHz", "float", "Uplink carrier frequency."),
    ("Uplink", ("eirp_gs_dbw",), "EIRP GS", "dBW", "float", "Ground station EIRP."),
    ("Uplink", ("gt_sat_dbk",), "G/T SAT", "dB/K", "float", "Satellite G/T."),
    ("Uplink", ("uplink_bitrate_mbps",), "Bit rate", "Mbps", "float", "Uplink coded bit rate."),
    ("Uplink", ("uplink_overhead",), "Overhead", "", "float", "Coded / information bit rate."),
    ("Uplink", ("uplink_rolloff",), "Roll-off factor", "", "float", ""),
    ("Uplink", ("uplink_spectral_efficiency_bpshz",), "Spectral efficiency", "b/s/Hz", "float", ""),
    ("Uplink", ("demod_loss_sat_db",), "Demodulator loss SAT", "dB", "float", ""),
    ("Uplink", ("uplink_other_attenuations_db",), "Other attenuations", "dB", "float", ""),
    ("Statistics", ("downlink_data_kbps", "instruments"), "Instrument data", "kbps", "float",
     "Payload data stored on board and sent on this downlink."),
    ("Statistics", ("downlink_data_kbps", "hktm"), "HKTM (stored)", "kbps", "float", "Housekeeping telemetry stored on board."),
    ("Statistics", ("downlink_data_kbps", "other"), "Other data", "kbps", "float", "E.g. events and system logs."),
    ("Statistics", ("downlink_data_kbps", "realtime"), "Real-time TM", "kbps", "float",
     "Sent live during the passes; takes part of the downlink bit rate."),
    ("Statistics", ("memory_gb",), "On-board memory", "GB", "float", "Memory of each satellite for these data."),
    ("Statistics", ("min_pass_min",), "Minimum pass duration", "min", "float", "Shorter contacts are not used."),
    ("Statistics", ("antenna_setup_min",), "Antenna setup", "min", "float", "Time to configure an antenna for another satellite."),
    ("Statistics", ("latency_at_pass_end",), "Latency at end of pass", "", "bool", "Data on ground when the dump is complete."),
    ("Statistics", ("product_length_min",), "Product length", "min", "float",
     "Ground products of this many minutes of on-board data (fixed blocks from 00:00 UTC); empty: none."),
    ("Statistics", ("product_latency_reference",), "Product latency from", "", "str",
     "data (each datum's sensing), start or end (of the product's sensing)."),
    ("Statistics", ("latency_requirement", "minutes"), "Latency requirement", "min", "float", ""),
    ("Statistics", ("latency_requirement", "percent"), "Latency requirement, share", "%", "float", ""),
    ("Statistics", ("statistics_stations",), "Ground stations", "", "list", "Network of the statistics."),
]

# Sections whose scalar values are listed with their keys (the static link budgets).
SECTIONS = {
    "fixed_elevation_downlink": "Static LB · Downlink",
    "fixed_elevation_uplink": "Static LB · Uplink",
}


def get_path(payload: dict, path: tuple[str, ...]) -> Any:
    node: Any = payload
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value)


def _kind_of(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, list):
        return "list"
    return "str"


def rows(payload: dict) -> list[dict]:
    """Every parameter as {group, path, label, unit, kind, description, value} (schema first, then the rest)."""
    result = []
    known = set()
    for group, path, label, unit, kind, description in SCHEMA:
        known.add(path)
        result.append(dict(group=group, path=path, label=label, unit=unit, kind=kind, description=description,
                           value=get_path(payload, path)))
    for section, group in SECTIONS.items():
        values = payload.get(section)
        if isinstance(values, dict):
            for key, value in values.items():
                if isinstance(value, dict):
                    continue  # e.g. tolerances: edited in the static link budget tab
                result.append(dict(group=group, path=(section, key), label=key, unit="", kind=_kind_of(value),
                                   description="", value=value))
    covered_tops = {p[0] for p in known} | set(SECTIONS)
    for key, value in payload.items():
        if key not in covered_tops and not isinstance(value, dict):
            result.append(dict(group="Other", path=(key,), label=key, unit="", kind=_kind_of(value),
                               description="", value=value))
    return result


def parse_value(text: str, kind: str) -> Any:
    """The value typed in ``text`` for a parameter of ``kind`` (None when empty). Raises ValueError."""
    text = text.strip()
    if not text:
        return None
    if kind == "float":
        return float(text.replace(",", "."))
    if kind == "int":
        return int(float(text))
    if kind == "bool":
        lowered = text.lower()
        if lowered in ("yes", "true", "1", "si", "sì", "on"):
            return True
        if lowered in ("no", "false", "0", "off"):
            return False
        raise ValueError(f"'{text}' is not yes or no")
    if kind == "list":
        return [part.strip() for part in text.split(",") if part.strip()]
    return text


def set_path(payload: dict, path: tuple[str, ...], value: Optional[Any]) -> None:
    """Set (or, with None, remove) the parameter at ``path``, creating the sections it needs."""
    node = payload
    for key in path[:-1]:
        if not isinstance(node.get(key), dict):
            if value is None:
                return
            node[key] = {}
        node = node[key]
    if value is None:
        node.pop(path[-1], None)
    else:
        node[path[-1]] = value


def standardize(payload: dict, custom_orbit: bool) -> dict:
    """The same parameters with every key of the schema, in the schema order.

    A key the mission does not use is written as null, so every
    parameters.json shows the same keys and where to fill them (e.g. the
    instrument data rate). The orbit section only for a custom-orbit mission.
    The static link budget sections and any other key are kept as they are,
    after the schema keys.
    """
    result: dict = {}
    for _group, path, _label, _unit, _kind, _description in SCHEMA:
        if path[0] == "orbit" and not custom_orbit and "orbit" not in payload:
            continue
        set_path_keep(result, path, get_path(payload, path))
    for key, value in payload.items():
        if key not in result:
            result[key] = value
        elif isinstance(value, dict) and isinstance(result[key], dict):
            for sub, sub_value in value.items():  # keys of a schema section not in the schema
                result[key].setdefault(sub, sub_value)
    return result


def set_path_keep(payload: dict, path: tuple[str, ...], value: Any) -> None:
    """Set the value at ``path`` (None included), creating the sections it needs."""
    node = payload
    for key in path[:-1]:
        if not isinstance(node.get(key), dict):
            node[key] = {}
        node = node[key]
    node[path[-1]] = value
