"""Download a satellite's TLE from the internet and check it.

Two public sources, both without login or API key, tried in this order:

* SatNOGS DB (db.satnogs.org): JSON API, TLE taken from Space-Track;
* CelesTrak (celestrak.org): the GP data in TLE format.

Only the standard library is used (urllib), so the system proxy settings
apply; on Windows these are the ones of the Internet options.
"""
from __future__ import annotations

import datetime
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
