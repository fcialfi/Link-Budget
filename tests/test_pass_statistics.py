import datetime

import numpy as np
import pytest

import pass_statistics as ps

UTC = datetime.timezone.utc


def _day(day, visible, margin=None, station=None):
    start = datetime.datetime(2026, 10, day, tzinfo=UTC)
    visible = np.array(visible, dtype=bool)
    return ps.DaySamples(
        start,
        np.where(visible, 20.0, -10.0),
        visible,
        None if margin is None else np.array(margin, float),
        None if station is None else np.array(station, dtype=object),
    )


def test_passes_capacity_and_gaps():
    # 6 samples of 10 minutes per "day" (step 600 s), two days.
    d1 = _day(1, [0, 1, 1, 0, 1, 0], [0, 4.0, 1.0, 0, 5.0, 0])
    d2 = _day(2, [1, 0, 0, 0, 0, 1], [3.0, 0, 0, 0, 0, 2.0])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6, min_margin_db=2.0)
    assert [p.contact_min for p in result.passes] == [20.0, 10.0, 10.0, 10.0]
    # Usable: margin >= 2 dB -> 1 of 2 samples in the first pass.
    assert result.passes[0].usable_min == 10.0
    assert result.passes[0].capacity_mb == pytest.approx(600 * 1e6 / 8e6)  # 75 MB
    assert result.passes[0].downloaded_mb == result.passes[0].capacity_mb  # no memory model
    assert result.passes[0].min_margin_db == 1.0
    day1, day2 = result.days
    assert day1.contacts == 2 and day2.contacts == 2
    assert day1.capacity_mb == pytest.approx(150.0)
    # Gaps: pass 1 ends at 00:30, pass 2 starts 00:40 -> 10 min; pass 2 ends 00:50, pass 3 starts next day 00:00.
    assert day1.max_gap_min == 10.0
    assert day2.max_gap_min == pytest.approx(24 * 60 - 50)
    assert result.summary["max_gap_min"] == pytest.approx(24 * 60 - 50)
    assert result.summary["downloaded_mb_total"] == pytest.approx(300.0)
    assert day1.on_board_mb is None


def test_pass_over_midnight_counts_for_its_first_day():
    d1 = _day(1, [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 3, 3])
    d2 = _day(2, [1, 0, 0, 0, 0, 0], [3, 0, 0, 0, 0, 0])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6)
    assert len(result.passes) == 1
    assert result.passes[0].contact_min == 30.0
    assert result.days[0].contacts == 1 and result.days[1].contacts == 0


def test_memory_model_limits_the_download_and_gives_the_latency():
    # Contact only in the last sample of each 6-sample day; the downlink could take 75 MB
    # per sample but only 6 samples x 600 s x 10 kbit/s = 4.5 MB are generated per day.
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    d2 = _day(2, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4)
    assert [d.downloaded_mb for d in result.days] == pytest.approx([4.5, 4.5])
    assert [d.capacity_mb for d in result.days] == pytest.approx([75.0, 75.0])
    assert [d.on_board_mb for d in result.days] == pytest.approx([0.0, 0.0])
    assert result.summary["generated_mb_per_day"] == pytest.approx(108.0)  # 10 kbit/s over 24 h
    # Oldest bit of the first dump: generated at 0, downloaded at the end of sample 6 -> 60 min.
    assert result.summary["max_latency_min"] == pytest.approx(60.0)
    # Data generated uniformly over the hour before the dump: mean latency 30 min.
    assert result.summary["mean_latency_min"] == pytest.approx(30.0)


def test_memory_fills_when_the_downlink_is_too_small():
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    # 1 Mbit/s generated: 450 MB per 6 samples, 75 MB downloaded.
    result = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e6)
    assert result.days[0].downloaded_mb == pytest.approx(75.0)
    assert result.days[0].on_board_mb == pytest.approx(375.0)
    assert result.summary["max_on_board_mb"] == pytest.approx(375.0)


def test_station_network():
    start = datetime.datetime(2026, 10, 1, tzinfo=UTC)
    a = ps.DaySamples(start, np.array([10.0, 30.0, -5.0, -5.0]), np.array([1, 1, 0, 0], bool), np.array([1.0, 6.0, np.nan, np.nan]))
    b = ps.DaySamples(start, np.array([-5.0, 20.0, 15.0, -5.0]), np.array([0, 1, 1, 0], bool), np.array([np.nan, 8.0, 4.0, np.nan]))
    combined = ps.combine_stations({"A": a, "B": b})
    assert list(combined.visible) == [True, True, True, False]
    assert list(combined.station) == ["A", "B", "B", ""]  # the best margin wins
    assert list(combined.margin_db[:3]) == [1.0, 8.0, 4.0]
    result = ps.aggregate([combined], step_s=60, info_bitrate_bps=1e6)
    assert len(result.passes) == 1
    assert result.passes[0].stations == "A + B"
    assert result.passes[0].contact_min == 3.0


def test_contact_analysis_only():
    result = ps.aggregate([_day(1, [0, 1, 0, 1, 1, 0])], step_s=600)
    assert result.days[0].downloaded_mb is None
    assert result.summary["passes"] == 2
    assert "downloaded_mb_total" not in result.summary


def test_compute_calls_each_day():
    seen = []

    def samples(day):
        seen.append(day)
        return ps.DaySamples(day, np.zeros(3), np.zeros(3, bool))

    ps.compute(samples, datetime.date(2026, 10, 30), 3, 600)
    assert [d.day for d in seen] == [30, 31, 1]
    with pytest.raises(ValueError):
        ps.compute(samples, datetime.date(2026, 10, 30), 0, 600)


def test_latency_requirement():
    # As in test_memory_model...: each dump takes the data of the hour before it.
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    result = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4,
                          latency_requirement=(45.0, 90.0))
    # Latencies uniform from 0 to 60 min: 75 % within 45 min, 90 % within 54 min.
    assert result.latency_fraction(45.0) == pytest.approx(0.75)
    assert result.summary["latency_within_pct"] == pytest.approx(75.0)
    assert result.summary["latency_at_required_pct_min"] == pytest.approx(54.0, abs=0.01)
    assert result.summary["latency_requirement_met"] is False
    assert result.latency_fraction(60.0) == pytest.approx(1.0)
