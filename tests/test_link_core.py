"""Regression tests of the pass analysis (link_core), on a frozen AWS case.

Reference: AWS over Svalbard on 14 October 2026, with the parameters and
antenna pattern of tests/data/aws_reference and the TLE below. The expected
values are those of the program at version 1.1.0 (the window gave the same
before the pass analysis moved to link_core).
"""
import datetime
import json
import os

import numpy as np
import pytest

import calculations
import link_core

DATA = os.path.join(os.path.dirname(__file__), "data", "aws_reference")
TLE = (
    "1 60543U 24149CD  26287.16919301  .00000688  00000-0  77235-4 0  9993",
    "2 60543  97.6685 354.1550 0003585  83.9816 276.1810 14.89818822117016",
)
SVALBARD = link_core.Station("Svalbard", 78.2232, 15.6267, 10.0, d_gs=3.0)
DAY = datetime.datetime(2026, 10, 14, tzinfo=datetime.timezone.utc)


@pytest.fixture(scope="module")
def params():
    with open(os.path.join(DATA, "parameters.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def aws_day(params):
    saved = (calculations.ANTENNA_PATTERN_ANGLES, calculations.ANTENNA_PATTERN_GAINS)
    calculations.load_antenna_pattern(os.path.join(DATA, "antenna_pattern.txt"))
    try:
        link = link_core.link_inputs_from_parameters(params)
        analysis = link_core.PassAnalysis(TLE, SVALBARD, link)
        frame = analysis.day(DAY)
        contacts = analysis.contacts(frame)
        yield link, analysis, frame, contacts
    finally:
        calculations.set_antenna_pattern(*saved)


def _at(frame, hh, mm, ss=0):
    t = DAY.replace(hour=hh, minute=mm, second=ss)
    return frame[frame["Time (UTC)"] == t].iloc[0]


def test_static_values_from_the_parameters(params):
    link = link_core.link_inputs_from_parameters(params)
    assert link.loss_kwargs["tx_chain_dbw"] == pytest.approx(4.5079, abs=1e-4)  # 3.3 W - 0.5 dB - VSWR 1.5
    assert link.loss_kwargs["link_losses_db"] == pytest.approx(2.6820, abs=1e-4)  # iono S4 + polarisation
    assert link.loss_kwargs["formatting_overhead"] == 1.00905
    assert link.info_bitrate_bps == pytest.approx(3.57e6 / 2.31493)
    assert (link.pfd_limit_low, link.pfd_limit_high) == (-154.0, -144.0)
    assert link.occupied_bandwidth_hz == pytest.approx(3.57e6)


def test_contacts(aws_day):
    _link, _analysis, _frame, contacts = aws_day
    assert len(contacts) == 15
    start, end, direction = contacts[0]
    assert (start.strftime("%H:%M:%S"), end.strftime("%H:%M:%S"), direction) == ("01:11:00", "01:15:00", "A")
    assert contacts[1][0].strftime("%H:%M:%S") == "02:50:30"


def test_link_budget_samples(aws_day):
    _link, _analysis, frame, _contacts = aws_day
    low = _at(frame, 1, 12)
    assert low["Elevation (°)"] == pytest.approx(6.71, abs=0.01)
    assert low["Eb/No (dB)"] == pytest.approx(13.01, abs=0.011)
    assert low["G/T GS (dB/K)"] == pytest.approx(9.79, abs=0.011)
    assert low["EIRP (dBW)"] == pytest.approx(5.01, abs=0.011)
    assert low["PFD Margin (dB)"] == pytest.approx(9.22, abs=0.011)
    high = _at(frame, 10, 58)
    assert high["Elevation (°)"] == pytest.approx(59.73, abs=0.01)
    assert high["Eb/No (dB)"] == pytest.approx(27.60, abs=0.011)
    assert high["EIRP (dBW)"] == pytest.approx(8.62, abs=0.011)
    assert high["PFD Margin (dB)"] == pytest.approx(4.72, abs=0.011)


def test_visibility_and_budget_rows(aws_day):
    _link, _analysis, frame, _contacts = aws_day
    visible = frame["Visible"].to_numpy()
    assert np.all(frame.loc[visible, "Elevation (°)"] >= 5.0)
    # The link budget is computed above the horizon only.
    assert frame.loc[frame["Elevation (°)"] <= 0, "Eb/No (dB)"].isna().all()
    assert frame.loc[visible, "Eb/No (dB)"].notna().all()


def test_horizon_mask_shortens_the_contacts(params):
    mask = (np.array([0.0, 90.0, 180.0, 270.0]), np.array([3.0, 12.0, 2.0, 4.0]))
    station = link_core.Station("Svalbard", 78.2232, 15.6267, 10.0, d_gs=3.0, mask=mask)
    plain = link_core.PassAnalysis(TLE, SVALBARD, None).day(DAY)
    masked = link_core.PassAnalysis(TLE, station, None).day(DAY)
    assert masked["Visible"].sum() < plain["Visible"].sum()
    assert masked["Min Elevation (°)"].max() == pytest.approx(12.0, abs=0.01)


def test_contact_analysis_without_link_budget():
    frame = link_core.PassAnalysis(TLE, SVALBARD, None).day(DAY)
    assert "Eb/No (dB)" not in frame
    assert len(link_core.PassAnalysis(TLE, SVALBARD, None).contacts(frame)) == 15


def test_command_line_lists_the_missions(capsys):
    import lbtool

    lbtool.main(["missions"])
    out = capsys.readouterr().out
    assert "AWS" in out and "EPS-Sterna" in out


def test_sun_noise_lowers_gt_and_ebno(params):
    # 17 May 2026 ~20:53:30 UTC: the Sun is about 0.2° from AWS seen from Svalbard.
    import datetime as dt

    link = link_core.link_inputs_from_parameters(params)
    quiet = link_core.link_inputs_from_parameters(dict(params, sun_noise=False))
    day = link_core.utc_midnight(dt.date(2026, 5, 17))
    with_sun = link_core.PassAnalysis(TLE, SVALBARD, link).day(day)
    without = link_core.PassAnalysis(TLE, SVALBARD, quiet).day(day)
    assert "Sun Noise (dB)" not in without and "Sun Separation (°)" in without
    sample = with_sun["Sun Noise (dB)"].idxmax()
    assert with_sun.at[sample, "Sun Separation (°)"] < 1.0
    assert with_sun.at[sample, "Sun Noise (dB)"] > 5.0
    drop = without.at[sample, "Eb/No (dB)"] - with_sun.at[sample, "Eb/No (dB)"]
    assert drop == pytest.approx(with_sun.at[sample, "Sun Noise (dB)"], abs=0.02)
    # Far from the Sun the noise is negligible.
    far = with_sun["Sun Separation (°)"] > 30
    assert with_sun.loc[far, "Sun Noise (dB)"].max() < 0.02


def test_scintillation_adds_no_sky_noise(params):
    # Tropospheric scintillation is a fade, not an absorption: it lowers C but not G/T.
    plain = link_core.link_inputs_from_parameters(dict(params, sun_noise=False))
    scint = link_core.link_inputs_from_parameters(dict(params, sun_noise=False))
    scint.scintillation = True
    a = link_core.PassAnalysis(TLE, SVALBARD, plain).day(DAY, budget="visible")
    b = link_core.PassAnalysis(TLE, SVALBARD, scint).day(DAY, budget="visible")
    seen = a["Visible"]
    assert (b.loc[seen, "Atmospheric Att (dB)"] > a.loc[seen, "Atmospheric Att (dB)"]).all()
    assert np.allclose(b.loc[seen, "G/T GS (dB/K)"], a.loc[seen, "G/T GS (dB/K)"], atol=0.011)
