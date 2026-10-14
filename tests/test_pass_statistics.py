import datetime

import numpy as np
import pytest

import pass_statistics as ps

UTC = datetime.timezone.utc


def _day(day, visible, margin=None):
    start = datetime.datetime(2026, 10, day, tzinfo=UTC)
    visible = np.array(visible, dtype=bool)
    return ps.DaySamples(start, np.where(visible, 20.0, -10.0), visible, None if margin is None else np.array(margin, float))


def test_passes_volume_and_gaps():
    # 6 samples of 10 minutes per "day" (step 600 s), two days.
    d1 = _day(1, [0, 1, 1, 0, 1, 0], [0, 4.0, 1.0, 0, 5.0, 0])
    d2 = _day(2, [1, 0, 0, 0, 0, 1], [3.0, 0, 0, 0, 0, 2.0])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6, min_margin_db=2.0, onboard_mb_per_day=100.0)
    assert [p.contact_min for p in result.passes] == [20.0, 10.0, 10.0, 10.0]
    # Usable: margin >= 2 dB -> 1 of 2 samples in the first pass.
    assert result.passes[0].usable_min == 10.0
    assert result.passes[0].volume_mb == pytest.approx(600 * 1e6 / 8e6)  # 75 MB
    assert result.passes[0].min_margin_db == 1.0
    day1, day2 = result.days
    assert day1.contacts == 2 and day2.contacts == 2
    assert day1.volume_mb == pytest.approx(150.0)
    # Gaps: pass 1 ends at 00:30, pass 2 starts 00:40 -> 10 min; pass 2 ends 00:50, pass 3 starts next day 00:00.
    assert day1.max_gap_min == 10.0
    assert day2.max_gap_min == pytest.approx(24 * 60 - 50)
    assert result.summary["max_gap_min"] == pytest.approx(24 * 60 - 50)
    # Backlog: 100 generated, 150 downloaded -> 0; day 2 downloads 150 again.
    assert day1.backlog_mb == 0.0 and day2.backlog_mb == 0.0
    assert result.summary["volume_mb_total"] == pytest.approx(300.0)


def test_pass_over_midnight_counts_for_its_first_day():
    d1 = _day(1, [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 3, 3])
    d2 = _day(2, [1, 0, 0, 0, 0, 0], [3, 0, 0, 0, 0, 0])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6)
    assert len(result.passes) == 1
    assert result.passes[0].contact_min == 30.0
    assert result.days[0].contacts == 1 and result.days[1].contacts == 0


def test_backlog_grows_when_the_downlink_is_too_small():
    d1 = _day(1, [1, 0, 0, 0, 0, 0], [5, 0, 0, 0, 0, 0])
    d2 = _day(2, [1, 0, 0, 0, 0, 0], [5, 0, 0, 0, 0, 0])
    result = ps.aggregate([d1, d2], step_s=600, info_bitrate_bps=1e6, onboard_mb_per_day=100.0)
    assert [d.backlog_mb for d in result.days] == [25.0, 50.0]
    assert result.summary["max_backlog_mb"] == 50.0


def test_contact_analysis_only():
    result = ps.aggregate([_day(1, [0, 1, 0, 1, 1, 0])], step_s=600)
    assert result.days[0].volume_mb is None
    assert result.summary["passes"] == 2
    assert "volume_mb_total" not in result.summary


def test_compute_calls_each_day():
    seen = []

    def samples(day):
        seen.append(day)
        return ps.DaySamples(day, np.zeros(3), np.zeros(3, bool))

    ps.compute(samples, datetime.date(2026, 10, 30), 3, 600)
    assert [d.day for d in seen] == [30, 31, 1]
    with pytest.raises(ValueError):
        ps.compute(samples, datetime.date(2026, 10, 30), 0, 600)
