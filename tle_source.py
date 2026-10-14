"""Download a satellite's TLE from the internet and check it.

Two public sources, both without login or API key, tried in this order:

* SatNOGS DB (db.satnogs.org): JSON API, TLE taken from Space-Track;
* CelesTrak (celestrak.org): the GP data in TLE format.

Only the standard library is used (urllib), so the system proxy settings
apply; on Windows these are the ones of the Internet options.
"""
from __future__ import annotations

import datetime
import math
import json
import urllib.request
from dataclasses import dataclass

SATNOGS_URL = "https://db.satnogs.org/api/tle/?norad_cat_id={norad_id}"
CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR={norad_id}&FORMAT=TLE"
TIMEOUT_S = 15
# A TLE older than this gives noticeably wrong pass times: warn the user.
MAX_AGE_DAYS = 7.0


@dataclass
class Tle:
    name: str
    line1: str
    line2: str
    source: str

    def text(self) -> str:
        """The TLE as written to tle.txt (name line, then the two lines)."""
        lines = [self.name] if self.name else []
        return "\n".join(lines + [self.line1, self.line2]) + "\n"


def checksum_ok(line: str) -> bool:
    """True when the last digit of a TLE line matches its modulo-10 checksum."""
    if len(line) < 69 or not line[68].isdigit():
        return False
    total = sum(int(c) if c.isdigit() else (1 if c == "-" else 0) for c in line[:68])
    return total % 10 == int(line[68])


def epoch(line1: str) -> datetime.datetime:
    """Epoch (UTC) of a TLE, from columns 19-32 of line 1 (YYDDD.DDDDDDDD)."""
    field = line1[18:32].strip()
    year = int(field[:2])
    year += 2000 if year < 57 else 1900
    day = float(field[2:])
    start = datetime.datetime(year, 1, 1, tzinfo=datetime.timezone.utc)
    return start + datetime.timedelta(days=day - 1)


def age_days(line1: str, now: datetime.datetime | None = None) -> float:
    """Days elapsed since the TLE epoch."""
    now = now or datetime.datetime.now(datetime.timezone.utc)
    return (now - epoch(line1)).total_seconds() / 86400.0


def _check(name: str, line1: str, line2: str, norad_id: int, source: str) -> Tle:
    line1, line2 = line1.strip(), line2.strip()
    if not (line1.startswith("1 ") and line2.startswith("2 ")):
        raise ValueError("the answer is not a TLE")
    if int(line1[2:7]) != norad_id:
        raise ValueError(f"TLE of NORAD {line1[2:7].strip()}, not {norad_id}")
    if not (checksum_ok(line1) and checksum_ok(line2)):
        raise ValueError("TLE checksum error")
    name = name.strip()
    if name.startswith("0 "):
        name = name[2:].strip()
    return Tle(name, line1, line2, source)


def _get(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Link-Budget-Tool"})
    with urllib.request.urlopen(request, timeout=TIMEOUT_S) as response:
        return response.read().decode("utf-8", errors="replace")


def _from_satnogs(norad_id: int) -> Tle:
    data = json.loads(_get(SATNOGS_URL.format(norad_id=norad_id)))
    if not data:
        raise ValueError(f"no TLE for NORAD {norad_id}")
    item = data[0]
    source = "SatNOGS"
    if item.get("tle_source"):
        source += f" ({item['tle_source']})"
    return _check(item.get("tle0", ""), item["tle1"], item["tle2"], norad_id, source)


def _from_celestrak(norad_id: int) -> Tle:
    lines = [line for line in _get(CELESTRAK_URL.format(norad_id=norad_id)).splitlines() if line.strip()]
    if len(lines) < 2 or not lines[-2].startswith("1 "):
        text = " ".join(lines)[:80] or "empty answer"
        raise ValueError(text)
    name = lines[-3] if len(lines) >= 3 else ""
    return _check(name, lines[-2], lines[-1], norad_id, "CelesTrak")


SOURCES = (("SatNOGS", _from_satnogs), ("CelesTrak", _from_celestrak))


def fetch_tle(norad_id: int) -> Tle:
    """Download the latest TLE of ``norad_id``, trying each source in turn.

    Raises ``RuntimeError`` listing why every source failed.
    """
    errors = []
    for name, source in SOURCES:
        try:
            return source(norad_id)
        except Exception as exc:  # network, HTTP, JSON or TLE errors
            errors.append(f"{name}: {exc or type(exc).__name__}")
    raise RuntimeError("Could not download the TLE:\n" + "\n".join(f"- {e}" for e in errors))


# ---------------------------------------------------------------------------
# Synthetic TLE for a mission not launched yet
# ---------------------------------------------------------------------------

MU_KM3_S2 = 398600.4418
EARTH_RADIUS_KM = 6378.137
J2 = 1.08262668e-3
# Rotation of the orbit plane needed for a sun-synchronous orbit (360° per tropical year).
SSO_NODE_RATE_RAD_S = 2 * math.pi / (365.2421897 * 86400.0)
SYNTHETIC_NORAD_ID = 99999


def sso_inclination_deg(altitude_km: float) -> float:
    """Inclination of the circular sun-synchronous orbit at ``altitude_km`` (J2 only)."""
    a = EARTH_RADIUS_KM + altitude_km
    n = math.sqrt(MU_KM3_S2 / a**3)
    cos_i = -SSO_NODE_RATE_RAD_S / (1.5 * n * J2 * (EARTH_RADIUS_KM / a) ** 2)
    if not -1.0 <= cos_i <= 1.0:
        raise ValueError(f"no sun-synchronous orbit at {altitude_km:g} km")
    return math.degrees(math.acos(cos_i))


def sun_right_ascension_deg(when: datetime.datetime) -> float:
    """Apparent right ascension of the Sun (low-precision formula, ~0.01°)."""
    jd = when.timestamp() / 86400.0 + 2440587.5
    n = jd - 2451545.0
    mean_lon = math.radians((280.460 + 0.9856474 * n) % 360.0)
    anomaly = math.radians((357.528 + 0.9856003 * n) % 360.0)
    ecl_lon = mean_lon + math.radians(1.915) * math.sin(anomaly) + math.radians(0.020) * math.sin(2 * anomaly)
    obliquity = math.radians(23.439 - 0.0000004 * n)
    ra = math.atan2(math.cos(obliquity) * math.sin(ecl_lon), math.cos(ecl_lon))
    return math.degrees(ra) % 360.0


def ltan_hours(text: str) -> float:
    """Parse a local time of the ascending node, "hh:mm" or decimal hours."""
    text = str(text).strip()
    if ":" in text:
        hours, minutes = text.split(":", 1)
        value = int(hours) + int(minutes) / 60.0
    else:
        value = float(text)
    if not 0.0 <= value < 24.0:
        raise ValueError(f"LTAN {text} is not between 00:00 and 24:00")
    return value


def raan_for_ltan_deg(ltan_h: float, when: datetime.datetime) -> float:
    """RAAN that puts the ascending node at local time ``ltan_h`` at ``when``."""
    return (sun_right_ascension_deg(when) + (ltan_h - 12.0) * 15.0) % 360.0


def ltan_of_tle(line1: str, line2: str) -> float:
    """Local time of the ascending node of a TLE, at its epoch (hours)."""
    raan = float(line2[17:25])
    return ((raan - sun_right_ascension_deg(epoch(line1))) / 15.0 + 12.0) % 24.0


def _with_checksum(line: str) -> str:
    total = sum(int(c) if c.isdigit() else (1 if c == "-" else 0) for c in line[:68])
    return line[:68] + str(total % 10)


def synthetic_tle(
    altitude_km: float,
    ltan_h: float,
    when: datetime.datetime,
    inclination_deg: float | None = None,
    name: str = "CUSTOM ORBIT",
    norad_id: int = SYNTHETIC_NORAD_ID,
) -> Tle:
    """TLE of a circular orbit, for the passes of a mission not launched yet.

    ``altitude_km`` is above the equatorial radius; without
    ``inclination_deg`` the orbit is sun-synchronous. The epoch is ``when``
    (use the start of the analysed day) and the satellite is at the
    ascending node at the epoch. No drag. The elements are taken as SGP4
    mean elements, which is accurate enough for preliminary pass geometry.
    """
    if not 100.0 <= altitude_km <= 40000.0:
        raise ValueError("altitude must be between 100 and 40000 km")
    if when.tzinfo is None:
        when = when.replace(tzinfo=datetime.timezone.utc)
    inclination = sso_inclination_deg(altitude_km) if inclination_deg is None else float(inclination_deg)
    if not 0.0 <= inclination <= 180.0:
        raise ValueError("inclination must be between 0 and 180°")
    a = EARTH_RADIUS_KM + altitude_km
    revs_per_day = math.sqrt(MU_KM3_S2 / a**3) * 86400.0 / (2 * math.pi)
    raan = raan_for_ltan_deg(ltan_h, when)

    start = datetime.datetime(when.year, 1, 1, tzinfo=datetime.timezone.utc)
    day = (when - start).total_seconds() / 86400.0 + 1.0
    line1 = (
        f"1 {norad_id:05d}U 00000A   {when.year % 100:02d}{day:012.8f}"
        f"  .00000000  00000-0  00000-0 0  999"
    )
    line2 = (
        f"2 {norad_id:05d} {inclination:8.4f} {raan:8.4f} 0001000   0.0000   0.0000 "
        f"{revs_per_day:11.8f}    1"
    )
    line1, line2 = _with_checksum(line1.ljust(68)), _with_checksum(line2.ljust(68))
    return Tle(name, line1, line2, "synthetic")
