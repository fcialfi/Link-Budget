import astropy.units as u
import pytest

import calculations
import static_report


@pytest.fixture
def budget(monkeypatch):
    def _fake(lat, lon, freq_ghz, elevation_deg, p, d_gs, hs, **kwargs):
        if kwargs.get("include_rain"):
            return (0.6, 0.2, 1.9, 0.3, 3.0)
        return (0.6, 0.0, 0.0, 0.0, 0.6)

    monkeypatch.setattr(calculations.itu, "atmospheric_attenuation_slant_path", _fake)
    return calculations.calculate_fixed_elevation_link_budget(
        freq=1.707 * u.GHz,
        elevation_deg=3.0,
        sat_altitude_km=628,
        lat_gs=78.2,
        lon_gs=15.6,
        alt_gs_km=0.01,
        d_gs=3.0,
        eirp=None,
        gt=7.0,
        demod_loss=1.0,
        bitrate=3570e3,
        overhead=2.29,
        other_att=0.0,
        pointing_loss_db=0.0,
        link_availability_pct=99.99,
        required_ebno=2.6,
        tx_power_w=3.3,
        antenna_circuit_loss_db=0.5,
        vswr=1.5,
        antenna_gain_dbi=0.0,
        ionospheric_loss_db=0.05,
        polarisation_loss_db=0.09,
        modulation_degradation_db=0.5,
        occupied_bandwidth_hz=5176.5e3,
        pfd_limit_dbw_m2_4khz=-154.0,
        formatting_overhead=1.02,
        tolerances=[{"name": "Tx Antenna Gain", "favourable_db": 1.5, "adverse_db": 1.5, "distribution": "TRI"}],
    )


EXTRA_INFO = {
    "station": "Svalbard (78.2232°, 15.6267°, 10 m)",
    "link_availability_pct": 99.99,
    "occupied_bandwidth_khz": 5176.5,
    "required_margin_db": 3.0,
    "bandwidth_missing": False,
}


def test_link_budget_rows_groups_and_checks(budget):
    rows = static_report.link_budget_rows(budget, required_margin_db=3.0)
    by_label = {row.label: row for row in rows}
    assert [row.group for row in rows][0] == "Transmitter"
    margin = by_label["Data Recovery Margin [dB]"]
    assert margin.clear == pytest.approx(budget["clear"]["margin_db"])
    assert margin.rain == pytest.approx(budget["rain_faded"]["margin_db"])
    assert margin.status(margin.clear) == "pass"
    assert margin.status(1.0) == "fail"
    assert by_label["MEAN – 3σ [dB]"].threshold == 0.0
    # Unchecked rows never get a status.
    assert by_label["EIRP [dBW]"].status(by_label["EIRP [dBW]"].clear) is None


def test_link_budget_rows_without_margin_requirement(budget):
    rows = static_report.link_budget_rows(budget)
    margin = next(row for row in rows if row.label == "Data Recovery Margin [dB]")
    assert margin.status(margin.clear) is None


def test_pfd_rows_and_notes(budget):
    rows = static_report.pfd_rows(budget, EXTRA_INFO)
    assert rows[-1].label == "PFD Margin [dB]"
    assert rows[-1].status() == "pass"
    assert static_report.pfd_notes(budget, EXTRA_INFO) == []
    assert static_report.pfd_notes(budget, {"bandwidth_missing": True})


def test_rich_text_uses_symbol_font_and_escapes():
    text = static_report._rich("MEAN – 3σ ≥ 0 & <x> ✓")
    assert '<font name="Symbol">σ</font>' in text
    assert '<font name="Symbol">≥</font>' in text
    assert "&amp;" in text and "&lt;x&gt;" in text
    assert "✓" not in text


def test_build_static_report_pdf(tmp_path, budget):
    pytest.importorskip("reportlab")
    path = tmp_path / "report.pdf"
    context = static_report.ReportContext(
        title="Fixed-Elevation Downlink Link Budget",
        direction="Downlink (S/C → GS)",
        station="Svalbard",
        elevation_deg=3.0,
        link_availability_pct=99.99,
        input_groups=[("Geometry & Requirements", [("Elevation Angle [deg]", "3")])],
    )
    static_report.build_static_report_pdf(str(path), budget, EXTRA_INFO, context)
    data = path.read_bytes()
    assert data.startswith(b"%PDF")
    assert b"/Title" in data
