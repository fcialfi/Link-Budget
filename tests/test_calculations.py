import pytest

np = pytest.importorskip('numpy')
skyfield = pytest.importorskip('skyfield.api')
itur = pytest.importorskip('itur')
from astropy import units as u
from skyfield.api import load, EarthSatellite, wgs84

from link_budget_tool.calculations import calculate_link_budget_parameters


def test_calculate_link_budget_basic():
    ts = load.timescale()
    # Simple TLE for ISS
    tle1 = "1 25544U 98067A   20029.54791435  .00000398  00000-0  11189-4 0  9999"
    tle2 = "2 25544  51.6435 231.2046 0007417  72.8980 338.7556 15.49536738211442"
    sat = EarthSatellite(tle1, tle2, "ISS", ts)
    gs = wgs84.latlon(0.0, 0.0, elevation_m=0)

    t = ts.utc(2020, 1, 1, 0, 0, 0)

    params = calculate_link_budget_parameters(
        t,
        sat,
        gs,
        2.2 * u.GHz,
        p=50.0,
        r001=0.1,
        d_gs=100.0,
        alt_gs=0.0,
        eirp_sat=30.0,
        gt_gs=20.0,
        demod_loss=0.0,
        bitrate=1e6,
        overhead=1.0,
        cisat_lin=None,
        other_att=0.0,
    )

    assert isinstance(params, dict)
    assert "Elevation (°)" in params
    assert "Visible" in params
