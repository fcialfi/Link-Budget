"""Unit tests for the pure calculation logic in ``calculations.py``.

These tests avoid the Tkinter GUI entirely and focus on the numerical core
of the link budget tool: ground station/antenna pattern file parsing,
Doppler shift, atmospheric attenuation error handling, and the link budget
formulas themselves.
"""

import os
from datetime import datetime, timezone

import astropy.units as u
import numpy as np
import pytest

import calculations


# --------------------------------------------------------------------------
# Ground station catalogue parsing
# --------------------------------------------------------------------------


def test_load_ground_stations_valid(tmp_path):
    gs_file = tmp_path / "ground_stations.txt"
    gs_file.write_text(
        "# comment line, should be ignored\n"
        "\n"
        "Station A, 10.0, 20.0, 100\n"
        "StationB,-5.5,30.25,0\n"
    )
    stations = calculations.load_ground_stations(str(gs_file))
    assert stations == {
        "Station A": (10.0, 20.0, 100.0),
        "StationB": (-5.5, 30.25, 0.0),
    }


def test_load_ground_stations_missing_file(tmp_path):
    missing = tmp_path / "nope.txt"
    with pytest.raises(FileNotFoundError):
        calculations.load_ground_stations(str(missing))


def test_load_ground_stations_wrong_field_count(tmp_path):
    gs_file = tmp_path / "ground_stations.txt"
    gs_file.write_text("OnlyThreeFields, 1.0, 2.0\n")
    with pytest.raises(ValueError, match="4 comma-separated fields"):
        calculations.load_ground_stations(str(gs_file))


def test_load_ground_stations_bad_numeric_value(tmp_path):
    gs_file = tmp_path / "ground_stations.txt"
    gs_file.write_text("Bad, not_a_number, 2.0, 3.0\n")
    with pytest.raises(ValueError, match="Invalid numeric value"):
        calculations.load_ground_stations(str(gs_file))


def test_load_ground_stations_empty_file_raises(tmp_path):
    gs_file = tmp_path / "ground_stations.txt"
    gs_file.write_text("# only comments\n\n")
    with pytest.raises(ValueError, match="No ground stations found"):
        calculations.load_ground_stations(str(gs_file))


def test_reload_ground_stations_updates_module_globals(tmp_path):
    gs_file = tmp_path / "ground_stations.txt"
    gs_file.write_text("Custom, 1.0, 2.0, 3\n")

    original_file = calculations.GROUND_STATIONS_FILE
    original_stations = calculations.GROUND_STATIONS
    try:
        stations = calculations.reload_ground_stations(str(gs_file))
        assert stations == {"Custom": (1.0, 2.0, 3.0)}
        assert calculations.GROUND_STATIONS == {"Custom": (1.0, 2.0, 3.0)}
        assert calculations.GROUND_STATIONS_FILE == os.path.abspath(str(gs_file))
        # load_ground_stations() with no argument must reflect the reload.
        assert calculations.load_ground_stations() == {"Custom": (1.0, 2.0, 3.0)}
    finally:
        calculations.GROUND_STATIONS_FILE = original_file
        calculations.GROUND_STATIONS = original_stations


# --------------------------------------------------------------------------
# Antenna pattern
# --------------------------------------------------------------------------


@pytest.fixture
def restore_antenna_pattern():
    angles = calculations.ANTENNA_PATTERN_ANGLES
    gains = calculations.ANTENNA_PATTERN_GAINS
    interp = calculations.PATTERN_INTERP
    max_gain = calculations.MAX_ANTENNA_GAIN
    yield
    calculations.ANTENNA_PATTERN_ANGLES = angles
    calculations.ANTENNA_PATTERN_GAINS = gains
    calculations.PATTERN_INTERP = interp
    calculations.MAX_ANTENNA_GAIN = max_gain


def test_antenna_pattern_zero_loss_at_boresight():
    assert calculations.antenna_pattern(0.0) == pytest.approx(0.0, abs=1e-9)


def test_antenna_pattern_scalar_returns_python_float():
    # Regression test: on NumPy 2.x, float() on a 1-element ndarray raises
    # TypeError instead of implicitly converting, which used to break every
    # scalar call to antenna_pattern().
    result = calculations.antenna_pattern(10.0)
    assert isinstance(result, float)
    assert not isinstance(result, np.ndarray)


def test_antenna_pattern_symmetric():
    assert calculations.antenna_pattern(30.0) == pytest.approx(
        calculations.antenna_pattern(-30.0)
    )


def test_antenna_pattern_clips_at_max_loss_outside_range():
    # Beyond the pattern's angle range, interpolation falls back to
    # fill_value=0 gain, so the pointing loss saturates at MAX_ANTENNA_GAIN.
    assert calculations.antenna_pattern(1000.0) == pytest.approx(
        calculations.MAX_ANTENNA_GAIN
    )


def test_antenna_pattern_array_input_preserves_shape():
    result = calculations.antenna_pattern(np.array([0.0, 65.0, 1000.0]))
    assert isinstance(result, np.ndarray)
    assert result.shape == (3,)
    assert result[0] == pytest.approx(0.0, abs=1e-9)


def test_load_antenna_pattern_replaces_globals(tmp_path, restore_antenna_pattern):
    pattern_file = tmp_path / "pattern.csv"
    pattern_file.write_text("-10,1.0\n0,2.0\n10,1.0\n")
    calculations.load_antenna_pattern(str(pattern_file))
    assert calculations.MAX_ANTENNA_GAIN == pytest.approx(2.0)
    assert calculations.antenna_pattern(0.0) == pytest.approx(0.0, abs=1e-9)


def test_load_antenna_pattern_rejects_single_column_file(tmp_path):
    pattern_file = tmp_path / "pattern.csv"
    pattern_file.write_text("1.0\n2.0\n3.0\n")
    with pytest.raises(ValueError):
        calculations.load_antenna_pattern(str(pattern_file))


# --------------------------------------------------------------------------
# Doppler shift
# --------------------------------------------------------------------------


class _FakeVector:
    def __init__(self, values):
        self.km = np.array(values, dtype=float)


class _FakeVelocity:
    def __init__(self, values):
        self.km_per_s = np.array(values, dtype=float)


class _FakeTopocentric:
    """Minimal stand-in for a Skyfield relative-position object."""

    def __init__(self, position_km, velocity_km_s):
        self.position = _FakeVector(position_km)
        self.velocity = _FakeVelocity(velocity_km_s)


def test_calculate_doppler_shift_approaching_is_positive():
    topo = _FakeTopocentric([1000.0, 0.0, 0.0], [-7.5, 0.0, 0.0])
    shift = calculations.calculate_doppler_shift(topo, 1.0 * u.GHz)
    assert shift > 0


def test_calculate_doppler_shift_receding_is_negative():
    topo = _FakeTopocentric([1000.0, 0.0, 0.0], [7.5, 0.0, 0.0])
    shift = calculations.calculate_doppler_shift(topo, 1.0 * u.GHz)
    assert shift < 0


def test_calculate_doppler_shift_matches_closed_form():
    topo = _FakeTopocentric([1000.0, 0.0, 0.0], [-7.5, 0.0, 0.0])
    freq = 2.0 * u.GHz
    shift_khz = calculations.calculate_doppler_shift(topo, freq)
    expected_khz = 7.5 / calculations.C_KM_S * freq.to(u.Hz).value / 1000.0
    assert shift_khz == pytest.approx(expected_khz)


def test_calculate_doppler_shift_purely_tangential_motion_is_zero():
    topo = _FakeTopocentric([1000.0, 0.0, 0.0], [0.0, 7.5, 0.0])
    shift = calculations.calculate_doppler_shift(topo, 1.0 * u.GHz)
    assert shift == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------
# Atmospheric attenuation
# --------------------------------------------------------------------------


def test_atmospheric_attenuation_falls_back_to_zero_on_error(monkeypatch):
    def _raise(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(calculations.itu, "atmospheric_attenuation_slant_path", _raise)
    result = calculations.atmospheric_attenuation(45.0, 8.0, 12.0, 0.1, 5.0, 0.1)
    assert result == 0.0


def test_atmospheric_attenuation_returns_nonnegative_float():
    result = calculations.atmospheric_attenuation(49.87, 8.65, 8.4, 0.1, 5.5, 0.144)
    assert isinstance(result, float)
    assert result >= 0.0


# --------------------------------------------------------------------------
# calculate_link_budget_parameters
# --------------------------------------------------------------------------


class _FakeTime:
    """Minimal stand-in for a Skyfield Time object."""

    def __init__(self, dt):
        self._dt = dt

    def utc_datetime(self):
        return self._dt


def _base_kwargs(**overrides):
    kwargs = dict(
        t_sky=_FakeTime(datetime(2024, 1, 1, tzinfo=timezone.utc)),
        sat=None,
        gs=None,
        freq=8.4 * u.GHz,
        p=0.1,
        d_gs=5.5,
        alt_gs=0.144,
        eirp_sat=20.0,
        gt_gs=15.0,
        demod_loss=1.0,
        bitrate=1e6,
        overhead=1.2,
        cisat_lin=None,
        other_att=0.5,
        atm_att=2.0,
        pre_alt_deg=45.0,
        pre_slant_range_km=1000.0,
        pre_doppler_khz=1.234,
        pre_off_boresight_angle=0.0,
    )
    kwargs.update(overrides)
    return kwargs


def test_link_budget_path_loss_matches_free_space_formula():
    params = calculations.calculate_link_budget_parameters(**_base_kwargs())
    expected_path_loss = 20 * np.log10(1000.0) + 20 * np.log10(8.4) + 92.45
    assert params["Path Loss (dB)"] == pytest.approx(expected_path_loss)


def test_link_budget_visible_above_min_elevation():
    params = calculations.calculate_link_budget_parameters(**_base_kwargs())
    assert params["Visible"] == "YES"


def test_link_budget_not_visible_below_min_elevation():
    params = calculations.calculate_link_budget_parameters(
        **_base_kwargs(pre_alt_deg=1.0)
    )
    assert params["Visible"] == "NO"


def test_link_budget_uses_precomputed_geometry_without_touching_sat_or_gs():
    # sat/gs are None here; if the function fell back to Skyfield calls it
    # would raise AttributeError, so a successful call proves the pre_*
    # shortcuts are honoured.
    params = calculations.calculate_link_budget_parameters(**_base_kwargs())
    assert params["Doppler Shift (kHz)"] == pytest.approx(1.234)
    assert params["Off Boresight Angle (°)"] == pytest.approx(0.0)


def test_link_budget_ebno_is_nan_when_bitrate_is_zero():
    params = calculations.calculate_link_budget_parameters(**_base_kwargs(bitrate=0))
    assert np.isnan(params["Eb/No (dB)"])


def test_link_budget_interference_reduces_cno():
    without_interference = calculations.calculate_link_budget_parameters(**_base_kwargs())
    with_interference = calculations.calculate_link_budget_parameters(
        **_base_kwargs(cisat_lin=10 ** (5.0 / 10.0))
    )
    assert (
        with_interference["C/(No+Io) (dBHz)"]
        < without_interference["C/(No+Io) (dBHz)"]
    )


def test_link_budget_uplink_fields_are_none_by_default():
    params = calculations.calculate_link_budget_parameters(**_base_kwargs())
    assert params["UL Path Loss (dB)"] is None
    assert params["UL Rx Power (dBW)"] is None
    assert params["UL Eb/No (dB)"] is None


def test_link_budget_uplink_computed_when_fully_specified():
    params = calculations.calculate_link_budget_parameters(
        **_base_kwargs(
            uplink_freq=6.0 * u.GHz,
            eirp_gs=50.0,
            gt_sat=-5.0,
            atm_att_ul=1.5,
        )
    )
    expected_ul_path_loss = 20 * np.log10(1000.0) + 20 * np.log10(6.0) + 92.45
    assert params["UL Path Loss (dB)"] == pytest.approx(expected_ul_path_loss)
    assert params["UL Rx Power (dBW)"] is not None
    assert params["UL Eb/No (dB)"] is not None
    assert np.isfinite(params["UL Eb/No (dB)"])


# --------------------------------------------------------------------------
# Integration: prepare_topocentric_data vs. calculate_link_budget_parameters
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def iss_context():
    from skyfield.api import load, EarthSatellite, wgs84

    tle1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9994"
    tle2 = "2 25544  51.6416 339.6194 0006317  58.1900 302.0968 15.49560328 32345"
    ts = load.timescale()
    sat = EarthSatellite(tle1, tle2, "ISS", ts)
    gs = wgs84.latlon(49.87, 8.65, elevation_m=144)
    return ts, sat, gs


def test_prepare_topocentric_data_shapes(iss_context):
    ts, sat, gs = iss_context
    times = ts.utc(2024, 1, 1, 0, 0, [0, 600, 1200])
    altitudes, slant_ranges, dopplers, off_boresight = calculations.prepare_topocentric_data(
        sat, gs, times, 8.4 * u.GHz
    )
    assert altitudes.shape == (3,)
    assert slant_ranges.shape == (3,)
    assert dopplers.shape == (3,)
    assert off_boresight.shape == (3,)
    assert np.all((off_boresight >= 0.0) & (off_boresight <= 180.0))


def test_link_budget_matches_prepare_topocentric_data(iss_context):
    """The two code paths the GUI relies on must agree with each other.

    ``run_analysis`` pre-computes geometry in bulk via
    ``prepare_topocentric_data`` and feeds the per-step values into
    ``calculate_link_budget_parameters``. Other call sites let the latter
    compute geometry itself. Both must produce the same numbers for the
    same time step.
    """

    ts, sat, gs = iss_context
    t_sky = ts.utc(2024, 1, 1, 12, 0, 0)
    freq = 8.4 * u.GHz

    params = calculations.calculate_link_budget_parameters(
        t_sky,
        sat,
        gs,
        freq,
        p=0.1,
        d_gs=5.5,
        alt_gs=0.144,
        eirp_sat=20.0,
        gt_gs=15.0,
        demod_loss=1.0,
        bitrate=1e6,
        overhead=1.2,
        cisat_lin=None,
        other_att=0.5,
    )

    altitudes, slant_ranges, dopplers, off_boresight = calculations.prepare_topocentric_data(
        sat, gs, ts.utc(2024, 1, 1, 12, 0, [0]), freq
    )

    assert params["Elevation (°)"] == pytest.approx(altitudes[0])
    assert params["Slant Range (km)"] == pytest.approx(slant_ranges[0])
    assert params["Doppler Shift (kHz)"] == pytest.approx(dopplers[0])
    assert params["Off Boresight Angle (°)"] == pytest.approx(off_boresight[0], abs=1e-6)


# --------------------------------------------------------------------------
# Fixed-elevation "static" link budget helpers
# --------------------------------------------------------------------------


def test_vswr_mismatch_loss_db_matches_reference_value():
    # EPS-STERNA-DDB-USERS reference sheet: VSWR 1.50:1 -> 0.18 dB.
    assert calculations.vswr_mismatch_loss_db(1.50) == pytest.approx(0.18, abs=1e-2)


def test_vswr_mismatch_loss_db_zero_for_perfect_match():
    assert calculations.vswr_mismatch_loss_db(1.0) == pytest.approx(0.0, abs=1e-9)


def test_power_flux_density_matches_reference_values():
    # EPS-STERNA-DDB-USERS reference sheet DESIGN column.
    result = calculations.power_flux_density(
        eirp_dbw=4.51, slant_range_km=2396.10, occupied_bandwidth_hz=5176.5e3
    )
    assert result["pfd_dbw_m2"] == pytest.approx(-134.07, abs=1e-2)
    assert result["pfd_dbw_m2_per_4khz"] == pytest.approx(-165.19, abs=1e-2)


def test_power_flux_density_omits_4khz_key_without_bandwidth():
    result = calculations.power_flux_density(eirp_dbw=4.51, slant_range_km=2396.10)
    assert "pfd_dbw_m2_per_4khz" not in result


@pytest.fixture
def _stub_atmospheric_contributions(monkeypatch):
    """Deterministic gas/cloud/rain/scintillation contributions for the fixed-budget tests."""

    def _fake(lat, lon, freq_ghz, elevation_deg, p, d_gs, hs, **kwargs):
        return (0.5, 0.0, 0.0, 0.0, 0.5)

    monkeypatch.setattr(calculations.itu, "atmospheric_attenuation_slant_path", _fake)


def test_fixed_elevation_link_budget_basic(_stub_atmospheric_contributions):
    budget = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.01,
        link_availability_pct=99.99,
        required_ebno=2.60,
    )

    assert budget["ionospheric_loss_db"] == 0.0
    assert budget["polarisation_loss_db"] == 0.0
    assert budget["multipath_loss_db"] == 0.0
    assert budget["modulation_degradation_db"] == 0.0
    assert "tx_power_dbw" not in budget
    assert "pfd_dbw_m2_per_4khz" not in budget
    assert budget["pfd_dbw_m2"] == pytest.approx(
        4.51 - 10 * np.log10(4 * np.pi) - 20 * np.log10(budget["slant_range_km"] * 1000.0)
    )
    clear = budget["clear"]
    assert clear["total_propagation_loss_db"] == pytest.approx(budget["path_loss_db"] + 0.5)
    assert clear["margin_db"] == pytest.approx(clear["ebno_db"] - 2.60)


def test_fixed_elevation_link_budget_extra_losses_reduce_rx_power(_stub_atmospheric_contributions):
    base = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
    )
    with_losses = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
        ionospheric_loss_db=0.01,
        polarisation_loss_db=0.04,
        multipath_loss_db=0.20,
        modulation_degradation_db=0.50,
    )

    total_extra_loss = 0.01 + 0.04 + 0.20 + 0.50
    assert with_losses["clear"]["rx_power_dbw"] == pytest.approx(
        base["clear"]["rx_power_dbw"] - 0.01 - 0.04 - 0.20
    )
    assert with_losses["clear"]["cno_dbhz"] == pytest.approx(
        base["clear"]["cno_dbhz"] - total_extra_loss
    )


def test_fixed_elevation_link_budget_tx_chain_breakdown(_stub_atmospheric_contributions):
    budget = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
        tx_power_w=3.30,
        antenna_circuit_loss_db=0.50,
        vswr=1.50,
        antenna_gain_dbi=0.0,
        occupied_bandwidth_hz=5176.5e3,
        pfd_limit_dbw_m2_4khz=-154.0,
    )

    assert budget["tx_power_dbw"] == pytest.approx(10 * np.log10(3.30))
    assert budget["vswr_loss_db"] == pytest.approx(0.18, abs=1e-2)
    assert budget["eirp_source"] == "tx_chain"
    assert budget["eirp_dbw"] == pytest.approx(
        budget["tx_power_dbw"] - 0.50 - budget["vswr_loss_db"]
    )
    assert budget["pfd_dbw_m2_per_4khz"] == pytest.approx(-165.19, abs=1e-2)
    assert budget["pfd_margin_db"] == pytest.approx(-154.0 - budget["pfd_dbw_m2_per_4khz"])


def test_fixed_elevation_link_budget_ebno_nan_without_bitrate(_stub_atmospheric_contributions):
    budget = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=0.0,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
    )
    assert np.isnan(budget["clear"]["ebno_db"])


def test_antenna_pointing_loss_db_matches_reference_value():
    # Reference sheet: 3 m dish, theta_3dB 3.97 deg, depointing 0.1 deg -> 0.01 dB.
    assert calculations.antenna_pointing_loss_db(0.1, 3.97) == pytest.approx(0.0076, abs=1e-4)
    assert calculations.antenna_pointing_loss_db(0.0, 3.97) == 0.0


def test_antenna_beamwidth_3db_deg_for_3m_dish_at_l_band():
    assert calculations.antenna_beamwidth_3db_deg(1.707, 3.0) == pytest.approx(4.10, abs=1e-2)


def test_polarisation_mismatch_loss_zero_for_perfect_circular_antennas():
    result = calculations.polarisation_mismatch_loss_db(0.0, 0.0)
    assert result["best"] == pytest.approx(0.0, abs=1e-12)
    assert result["worst"] == pytest.approx(0.0, abs=1e-12)


def test_polarisation_mismatch_loss_worst_exceeds_best():
    result = calculations.polarisation_mismatch_loss_db(3.0, 1.0)
    assert 0.0 < result["best"] < result["worst"]


@pytest.mark.parametrize(
    "dist, mean, variance",
    [
        # Reference sheet: antenna gain 0 dBi +/-1.5 dB TRI -> variance 0.38.
        ("TRI", 0.0, 0.375),
        # Reference sheet: effective gain range 3 dB UNI -> variance 0.75.
        ("UNI", 0.0, 0.75),
        ("GAU", 0.0, 0.25),
    ],
)
def test_tolerance_statistics_symmetric(dist, mean, variance):
    stats = calculations.tolerance_statistics(1.5, 1.5, dist)
    assert stats["mean_db"] == pytest.approx(mean)
    assert stats["variance_db2"] == pytest.approx(variance)


def test_tolerance_statistics_asymmetric_triangular_mean_shift():
    stats = calculations.tolerance_statistics(0.0, 1.5, "TRI")
    assert stats["mean_db"] == pytest.approx(-0.5)
    assert stats["variance_db2"] == pytest.approx(1.5**2 / 18)


def test_tolerance_statistics_rejects_unknown_distribution():
    with pytest.raises(ValueError):
        calculations.tolerance_statistics(1.0, 1.0, "FOO")


def test_margin_statistics_combines_contributions():
    stats = calculations.margin_statistics(
        5.0,
        [
            {"name": "gain", "favourable_db": 1.5, "adverse_db": 1.5, "distribution": "TRI"},
            {"name": "atm", "favourable_db": 0.0, "adverse_db": 2.0, "distribution": "UNI"},
        ],
    )
    assert stats["favourable_db"] == pytest.approx(6.5)
    assert stats["adverse_db"] == pytest.approx(1.5)
    assert stats["mean_db"] == pytest.approx(5.0 - 1.0)
    assert stats["variance_db2"] == pytest.approx(0.375 + 4.0 / 12)
    assert stats["mean_minus_3sigma_db"] == pytest.approx(4.0 - 3 * np.sqrt(0.375 + 4.0 / 12))
    assert stats["worst_case_rss_db"] == pytest.approx(5.0 - np.sqrt(1.5**2 + 2.0**2))
    assert [c["name"] for c in stats["contributions"]] == ["gain", "atm"]


def test_fixed_elevation_link_budget_reproduces_reference_sheet(monkeypatch):
    """Nominal column of the reference S/C -> E/S static link budget sheet."""

    monkeypatch.setattr(
        calculations.itu,
        "atmospheric_attenuation_slant_path",
        lambda *args, **kwargs: (3.0, 0.0, 0.0, 0.0, 3.0),
    )
    budget = calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=3.0,
        sat_altitude_km=627.99,
        lat_gs=78.2,
        lon_gs=15.4,
        alt_gs_km=0.0,
        d_gs=3.0,
        eirp=4.51,
        gt=7.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=3570 / 1559,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
        required_ebno=2.60,
        tx_power_w=3.30,
        antenna_circuit_loss_db=0.50,
        vswr=1.50,
        antenna_gain_dbi=0.0,
        ionospheric_loss_db=0.05,
        polarisation_loss_db=0.09,
        modulation_degradation_db=0.50,
        occupied_bandwidth_hz=5176.5e3,
        pfd_limit_dbw_m2_4khz=-154.0,
        formatting_overhead=1590.13 / 1559,
        tolerances=[{"name": "gain", "favourable_db": 1.5, "adverse_db": 1.5, "distribution": "TRI"}],
    )

    assert budget["slant_range_km"] == pytest.approx(2584.52, abs=0.1)
    assert budget["path_loss_db"] == pytest.approx(165.34, abs=1e-2)
    assert budget["effective_gain_dbi"] == pytest.approx(-0.68, abs=1e-2)
    assert budget["pfd_dbm_m2"] == pytest.approx(-104.73, abs=1e-2)
    assert budget["pfd_dbw_m2_per_4khz"] == pytest.approx(-165.85, abs=1e-2)
    assert budget["pfd_margin_db"] == pytest.approx(11.85, abs=1e-2)
    assert budget["bitrate_dbhz"] == pytest.approx(62.01, abs=1e-2)
    clear = budget["clear"]
    assert clear["total_propagation_loss_db"] == pytest.approx(168.48, abs=1e-2)
    assert clear["received_cno_dbhz"] == pytest.approx(71.63, abs=2e-2)
    assert clear["ebno_db"] == pytest.approx(8.12, abs=2e-2)
    assert clear["margin_db"] == pytest.approx(5.52, abs=2e-2)
    stats = clear["margin_statistics"]
    assert stats["favourable_db"] == pytest.approx(7.01, abs=2e-2)
    assert stats["adverse_db"] == pytest.approx(4.01, abs=2e-2)


def test_fixed_elevation_link_budget_polarisation_and_rx_pointing(_stub_atmospheric_contributions):
    common = dict(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
    )
    base = calculations.calculate_fixed_elevation_link_budget(**common)
    derived = calculations.calculate_fixed_elevation_link_budget(
        **common,
        polarisation_loss_db=None,
        tx_axial_ratio_db=3.0,
        rx_axial_ratio_db=1.0,
        rx_depointing_deg=0.5,
    )

    worst = calculations.polarisation_mismatch_loss_db(3.0, 1.0)["worst"]
    assert derived["polarisation_loss_db"] == pytest.approx(worst)
    assert derived["rx_beamwidth_3db_deg"] == pytest.approx(
        calculations.antenna_beamwidth_3db_deg(1.707, 3.0)
    )
    pointing = calculations.antenna_pointing_loss_db(0.5, derived["rx_beamwidth_3db_deg"])
    assert derived["rx_pointing_loss_db"] == pytest.approx(pointing)
    assert derived["clear"]["rx_power_dbw"] == pytest.approx(
        base["clear"]["rx_power_dbw"] - worst - pointing
    )


def test_fixed_elevation_link_budget_formatting_overhead_lowers_ebno(_stub_atmospheric_contributions):
    common = dict(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        eirp=4.51,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
    )
    base = calculations.calculate_fixed_elevation_link_budget(**common)
    formatted = calculations.calculate_fixed_elevation_link_budget(**common, formatting_overhead=1.02)
    assert formatted["clear"]["ebno_db"] == pytest.approx(base["clear"]["ebno_db"] - 10 * np.log10(1.02))
    assert "margin_statistics" not in formatted["clear"]


def test_fixed_elevation_link_budget_tx_chain_eirp_drives_rx_power(_stub_atmospheric_contributions):
    common = dict(
        freq=1.707 * u.GHz,
        elevation_deg=5.0,
        sat_altitude_km=628,
        lat_gs=78.9,
        lon_gs=11.9,
        alt_gs_km=0.1,
        d_gs=3.0,
        gt=5.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
    )
    from_input = calculations.calculate_fixed_elevation_link_budget(eirp=4.51, **common)
    from_chain = calculations.calculate_fixed_elevation_link_budget(
        eirp=None, tx_power_w=3.30, antenna_circuit_loss_db=0.50, vswr=1.50, antenna_gain_dbi=0.0, **common
    )
    assert from_input["eirp_source"] == "input"
    assert from_chain["eirp_source"] == "tx_chain"
    assert from_chain["clear"]["rx_power_dbw"] == pytest.approx(
        from_input["clear"]["rx_power_dbw"] + from_chain["eirp_dbw"] - 4.51
    )
    with pytest.raises(ValueError):
        calculations.calculate_fixed_elevation_link_budget(eirp=None, **common)
    with pytest.raises(ValueError):
        calculations.calculate_fixed_elevation_link_budget(eirp=None, tx_power_w=3.3, **common)
