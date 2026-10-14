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


def test_all_sources_fail(monkeypatch):
    def fake_get(url):
        if "satnogs" in url:
            return "[]"
        raise OSError("timed out")

    monkeypatch.setattr(tle_source, "_get", fake_get)
    with pytest.raises(RuntimeError) as err:
        tle_source.fetch_tle(60543)
    assert "SatNOGS: no TLE" in str(err.value)
    assert "N2YO: timed out" in str(err.value)


def test_wrong_satellite_rejected(monkeypatch):
    answer = [{"tle0": "X", "tle1": L1, "tle2": L2}]
    monkeypatch.setattr(tle_source, "_get", lambda url: json.dumps(answer) if "satnogs" in url else "")
    with pytest.raises(RuntimeError, match="not 12345"):
        tle_source.fetch_tle(12345)


def test_sso_inclination():
    assert tle_source.sso_inclination_deg(628) == pytest.approx(97.90, abs=0.02)
    assert tle_source.sso_inclination_deg(800) == pytest.approx(98.60, abs=0.02)


def test_ltan_parsing():
    assert tle_source.ltan_hours("22:20") == pytest.approx(22 + 1 / 3)
    assert tle_source.ltan_hours("9.5") == 9.5
    with pytest.raises(ValueError):
        tle_source.ltan_hours("25:00")


def test_aws_ltan_from_its_tle():
    # The ascending node of AWS is at about 22:20 local time.
    assert tle_source.ltan_of_tle(L1, L2) == pytest.approx(22.33, abs=0.05)


def test_synthetic_tle_is_valid_and_has_the_requested_orbit():
    when = datetime.datetime(2026, 10, 14, tzinfo=datetime.timezone.utc)
    tle = tle_source.synthetic_tle(628, tle_source.ltan_hours("22:20"), when, name="X")
    assert len(tle.line1) == len(tle.line2) == 69
    assert tle_source.checksum_ok(tle.line1) and tle_source.checksum_ok(tle.line2)
    assert tle_source.epoch(tle.line1) == when
    assert float(tle.line2[8:16]) == pytest.approx(97.8983, abs=1e-4)
    assert tle_source.ltan_of_tle(tle.line1, tle.line2) == pytest.approx(22 + 1 / 3, abs=1e-3)


def test_nodal_period_of_the_eps_sterna_reference_orbit():
    # EPS-Sterna reference orbit: a = 6973.702 km, e = 0.001214, i = 97.794°, nodal period 96.716 min.
    assert tle_source.nodal_period_min(6973.702, 0.001214, 97.794) == pytest.approx(96.716, abs=0.002)


def test_synthetic_tle_flies_the_eps_sterna_orbit():
    skyfield = pytest.importorskip("skyfield.api")
    np = pytest.importorskip("numpy")
    when = datetime.datetime(2026, 10, 14, tzinfo=datetime.timezone.utc)
    tle = tle_source.synthetic_tle(
        595.6, tle_source.ltan_hours("15:30"), when, 97.794, eccentricity=0.001214, arg_perigee_deg=90
    )
    assert tle.line2[26:33] == "0012140"
    ts = skyfield.load.timescale()
    sat = skyfield.EarthSatellite(tle.line1, tle.line2, "X", ts)
    times = ts.utc(2026, 10, 14, 0, np.arange(0, 1440, 0.25))
    position = sat.at(times)
    # Mean altitude over the ellipsoid 610.5 km (table of the reference orbit).
    assert skyfield.wgs84.height_of(position).km.mean() == pytest.approx(610.5, abs=1.0)
    # Nodal period from the ascending node crossings: 96.716 min, i.e. the passes do not drift.
    z = position.position.km[2]
    minutes = np.arange(0, 1440, 0.25)
    k = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]
    crossings = minutes[k] + 0.25 * (-z[k]) / (z[k + 1] - z[k])
    assert np.diff(crossings).mean() == pytest.approx(96.716, abs=0.01)
    # At the epoch the satellite is at the ascending node (mean elements: within ~0.5°).
    assert abs(z[0]) < 60


def test_synthetic_tle_with_inclination_and_bad_inputs():
    when = datetime.datetime(2026, 1, 1, tzinfo=datetime.timezone.utc)
    tle = tle_source.synthetic_tle(800, 10.0, when, inclination_deg=55.0)
    assert float(tle.line2[8:16]) == 55.0
    with pytest.raises(ValueError):
        tle_source.synthetic_tle(50, 10.0, when)


def test_n2yo_page_when_satnogs_does_not_have_the_satellite(monkeypatch):
    # A page with a real TLE of AWS, shown the way N2YO does (lines separated by <br>).
    page = f"<html><head><title>ARCTIC WEATHER SATELLITE</title></head><body><pre>{L1}<br/>{L2}</pre></body></html>"
    monkeypatch.setattr(tle_source, "_get", lambda url: "[]" if "satnogs" in url else page)
    tle = tle_source.fetch_tle(60543)
    assert tle.source == "N2YO"
    assert (tle.line1, tle.line2) == (L1, L2)
    assert tle.name == "ARCTIC WEATHER SATELLITE"


def test_n2yo_page_without_tle(monkeypatch):
    monkeypatch.setattr(tle_source, "_get", lambda url: "[]" if "satnogs" in url else "<html>nothing</html>")
    with pytest.raises(RuntimeError, match="N2YO: no TLE of NORAD 60543 on the page"):
        tle_source.fetch_tle(60543)
