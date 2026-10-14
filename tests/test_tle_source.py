import datetime
import json

import pytest

import tle_source

L1 = "1 60543U 24149CD  26287.16919301  .00000688  00000-0  77235-4 0  9993"
L2 = "2 60543  97.6685 354.1550 0003585  83.9816 276.1810 14.89818822117016"


def test_checksum():
    assert tle_source.checksum_ok(L1)
    assert tle_source.checksum_ok(L2)
    assert not tle_source.checksum_ok(L1[:-1] + "0")


def test_epoch_and_age():
    epoch = tle_source.epoch(L1)
    assert (epoch.year, epoch.month, epoch.day, epoch.hour) == (2026, 10, 14, 4)
    later = epoch + datetime.timedelta(days=2)
    assert tle_source.age_days(L1, now=later) == pytest.approx(2.0)


def test_satnogs_answer(monkeypatch):
    answer = [{"tle0": "0 ARCTIC WEATHER SATELLITE", "tle1": L1, "tle2": L2, "tle_source": "Space-Track.org"}]
    monkeypatch.setattr(tle_source, "_get", lambda url: json.dumps(answer))
    tle = tle_source.fetch_tle(60543)
    assert tle.name == "ARCTIC WEATHER SATELLITE"
    assert (tle.line1, tle.line2) == (L1, L2)
    assert tle.source == "SatNOGS (Space-Track.org)"
    assert tle.text() == f"ARCTIC WEATHER SATELLITE\n{L1}\n{L2}\n"


def test_falls_back_to_celestrak(monkeypatch):
    def fake_get(url):
        if "satnogs" in url:
            raise OSError("timed out")
        return f"AWS\r\n{L1}\r\n{L2}\r\n"

    monkeypatch.setattr(tle_source, "_get", fake_get)
    tle = tle_source.fetch_tle(60543)
    assert tle.source == "CelesTrak"
    assert tle.name == "AWS"


def test_all_sources_fail(monkeypatch):
    def fake_get(url):
        if "satnogs" in url:
            return "[]"
        return "No GP data found"

    monkeypatch.setattr(tle_source, "_get", fake_get)
    with pytest.raises(RuntimeError) as err:
        tle_source.fetch_tle(60543)
    assert "SatNOGS: no TLE" in str(err.value)
    assert "CelesTrak: No GP data found" in str(err.value)


def test_wrong_satellite_rejected(monkeypatch):
    answer = [{"tle0": "X", "tle1": L1, "tle2": L2}]
    monkeypatch.setattr(tle_source, "_get", lambda url: json.dumps(answer) if "satnogs" in url else "")
    with pytest.raises(RuntimeError, match="not 12345"):
        tle_source.fetch_tle(12345)
