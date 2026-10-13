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
        tolerances=[{"name": "Tx Antenna Gain", "key": "tx_gain", "favourable_db": 1.5, "adverse_db": 1.5, "distribution": "TRI"}],
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


def test_analysis_notes():
    assert static_report.analysis_notes({"elevation_deg": 5.0}) == []
    low = static_report.analysis_notes({"elevation_deg": 3.0})
    assert len(low) == 1 and "below 5°" in low[0]
    ignored = static_report.analysis_notes({"elevation_deg": 10.0, "tolerances_ignored": True})
    assert len(ignored) == 1 and "Required Eb/No" in ignored[0]


@pytest.mark.parametrize("applicable", [True, False])
def test_report_pfd_section_only_when_applicable(tmp_path, monkeypatch, budget, applicable):
    rl_config = pytest.importorskip("reportlab.rl_config")
    # Uncompressed content streams, so the page text can be searched directly.
    monkeypatch.setattr(rl_config, "pageCompression", 0)
    path = tmp_path / "report.pdf"
    context = static_report.ReportContext(
        title="Fixed-Elevation Link Budget",
        direction="Downlink (S/C → GS)" if applicable else "Uplink (GS → S/C)",
        station="Svalbard",
        elevation_deg=3.0,
        link_availability_pct=99.99,
    )
    extra = dict(EXTRA_INFO, pfd_applicable=applicable, elevation_deg=3.0)
    static_report.build_static_report_pdf(str(path), budget, extra, context)
    data = path.read_bytes()
    assert (b"Power Flux Density" in data) is applicable
    assert (b"PFD Margin" in data) is applicable
    assert b"below 5" in data


def test_other_attenuations_row_only_when_set(budget):
    labels = [row.label for row in static_report.link_budget_rows(budget)]
    assert "Other Attenuations [dB]" not in labels
    with_other = dict(budget, other_attenuation_db=0.5)
    rows = {row.label: row for row in static_report.link_budget_rows(with_other)}
    assert rows["Other Attenuations [dB]"].clear == 0.5


def test_analysis_notes_for_scintillation_and_other_attenuations():
    budget = {
        "scintillation_included": True,
        "tropospheric_scintillation_valid": False,
        "frequency_mhz": 1707.0,
        "ionospheric_scintillation": {"valid": False},
        "other_attenuation_db": 3.0,
    }
    notes = static_report.analysis_notes({"elevation_deg": 10.0}, budget)
    assert any("Tropospheric scintillation" in n for n in notes)
    assert any("S4 ≥ 0.4" in n for n in notes)
    assert any("counts the same effects twice" in n for n in notes)
    clean = dict(budget, scintillation_included=False, ionospheric_scintillation={"valid": True}, other_attenuation_db=0.5)
    assert static_report.analysis_notes({"elevation_deg": 10.0}, clean) == []


def test_info_rows_show_zone_and_s4(budget):
    with_scint = dict(
        budget,
        ionospheric_scintillation=calculations.ionospheric_scintillation_fade_db(0.3, 1707.0),
    )
    labels = dict(static_report.info_rows(with_scint, {}))
    assert "polar cap" in labels["Geomagnetic Latitude / Ionospheric Zone"]
    assert labels["Ionospheric S4 at Link Frequency"].startswith("0.266")


def test_tolerance_rows_propagate_to_the_margin(budget):
    rows = {row.label: row for row in static_report.tolerance_rows(budget, "rain_faded")}
    rain = budget["rain_faded"]
    # Only the Tx antenna gain has a tolerance (+/-1.5 dB TRI) in the fixture.
    eirp = rows["EIRP [dBW]"]
    assert eirp.nominal == pytest.approx(budget["eirp_dbw"])
    assert (eirp.favourable, eirp.adverse) == pytest.approx((budget["eirp_dbw"] + 1.5, budget["eirp_dbw"] - 1.5))
    assert eirp.variance == pytest.approx(0.375) and eirp.distribution == "TRI"
    margin = rows["Data Recovery Margin [dB]"]
    assert margin.nominal == pytest.approx(rain["margin_db"])
    assert margin.adverse == pytest.approx(rain["margin_statistics"]["adverse_db"])
    assert "Total Propagation Loss [dB]" not in rows  # no propagation tolerance


def test_tolerance_rows_losses_move_the_other_way(budget):
    stats = calculations.margin_statistics(
        budget["rain_faded"]["margin_db"],
        [{"name": "Prop", "key": "propagation", "favourable_db": 0.5, "adverse_db": 1.0, "distribution": "UNI"}],
    )
    with_prop = dict(budget, rain_faded=dict(budget["rain_faded"], margin_statistics=stats))
    loss = {row.label: row for row in static_report.tolerance_rows(with_prop)}["Total Propagation Loss [dB]"]
    nominal = budget["rain_faded"]["total_propagation_loss_db"]
    assert (loss.favourable, loss.adverse) == pytest.approx((nominal - 0.5, nominal + 1.0))
    assert loss.mean == pytest.approx(nominal + 0.25)


def test_summary_checks_and_verdict(budget):
    checks = static_report.summary_checks(budget, EXTRA_INFO)
    labels = [check.label for check in checks]
    assert labels[0] == "Received Eb/No [dB]"
    assert labels[-1] == "PFD Margin [dB]"
    margin = next(check for check in checks if check.label == "Data Recovery Margin [dB]")
    assert margin.rain == pytest.approx(budget["rain_faded"]["margin_db"])
    assert margin.requirement == "≥ 3 dB"
    assert static_report.summary_verdict(checks)[0] == "pass"

    # A margin requirement above the rain-faded margin fails the whole check.
    strict = dict(EXTRA_INFO, required_margin_db=budget["rain_faded"]["margin_db"] + 0.1)
    overall, text = static_report.summary_verdict(static_report.summary_checks(budget, strict))
    assert overall == "fail"
    assert text == "Requirement not met: Data Recovery Margin."

    # No PFD check on uplinks.
    uplink = dict(EXTRA_INFO, pfd_applicable=False)
    assert "PFD Margin [dB]" not in [c.label for c in static_report.summary_checks(budget, uplink)]


def test_summary_verdict_without_requirements():
    assert static_report.summary_verdict([])[0] is None
