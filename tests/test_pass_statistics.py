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


def _sat_day(visible, margin=5.0):
    start = datetime.datetime(2026, 10, 1, tzinfo=UTC)
    visible = np.array(visible, dtype=bool)
    return ps.DaySamples(start, np.where(visible, 20.0, -5.0), visible, np.where(visible, margin, np.nan))


def test_constellation_station_conflict():
    # Two satellites seen together by the only station in samples 1-2: one waits.
    day = {
        "S1": {"G": _sat_day([0, 1, 1, 0, 0, 0])},
        "S2": {"G": _sat_day([0, 1, 1, 1, 0, 0])},
    }
    result = ps.simulate_constellation([day], step_s=60, info_bitrate_bps=1e6, generated_bps=1e4)
    s1, s2 = result.satellites
    # Same memory: S1 (first) gets the station in sample 1 and keeps it while its pass lasts
    # (sample 2), although S2 then has more data; S2 is served in sample 3 only.
    assert s1["served_min_per_day"] == 2.0 and s2["served_min_per_day"] == 1.0
    assert s2["conflict_min_per_day"] == 2.0
    assert s1["conflict_min_per_day"] + s2["conflict_min_per_day"] == 2.0
    assert result.stations[0]["served_min_per_day"] == 3.0
    assert result.stations[0]["contention_min_per_day"] == 2.0
    assert result.summary["satellites"] == 2


def test_constellation_two_stations_serve_both():
    day = {
        "S1": {"A": _sat_day([0, 1, 1, 0]), "B": _sat_day([0, 1, 1, 0], margin=3.0)},
        "S2": {"A": _sat_day([0, 1, 1, 0], margin=6.0), "B": _sat_day([0, 0, 0, 0])},
    }
    result = ps.simulate_constellation([day], step_s=60, info_bitrate_bps=1e6, generated_bps=1e4)
    assert result.summary["conflict_min_per_day"] == 0.0
    by_station = {s["station"]: s for s in result.stations}
    assert by_station["A"]["served_min_per_day"] == 2.0
    assert by_station["B"]["served_min_per_day"] == 2.0
    # Data generated before each dump: everything generated up to the last dump is on ground.
    assert result.summary["downloaded_mb_per_day"] == pytest.approx(2 * 3 * 60 * 1e4 / 8e6)


def test_constellation_station_with_two_antennas_serves_both():
    day = {
        "S1": {"G": _sat_day([0, 1, 1, 0, 0, 0])},
        "S2": {"G": _sat_day([0, 1, 1, 1, 0, 0])},
    }
    result = ps.simulate_constellation([day], step_s=60, info_bitrate_bps=1e6, generated_bps=1e4,
                                       antennas={"G": 2})
    assert result.summary["conflict_min_per_day"] == 0.0
    station = result.stations[0]
    assert station["antennas"] == 2
    assert station["served_min_per_day"] == 5.0
    assert station["contention_min_per_day"] == 0.0
    assert station["busy_pct"] == pytest.approx(5 / 12 * 100)


def test_on_board_series_follows_the_memory():
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    result = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4)
    times, on_board = result.on_board_series
    assert len(times) == 6 and times[0] == d1.start
    # 0.75 MB generated per sample, all of it downloaded in the last one.
    assert on_board == pytest.approx([0.75, 1.5, 2.25, 3.0, 3.75, 0.0])
    assert ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6).on_board_series is None


def test_constellation_antenna_setup_time():
    # S1 in view in samples 1-2, S2 in samples 3-6; 2 min to reconfigure the antenna.
    day = {
        "S1": {"G": _sat_day([0, 1, 1, 0, 0, 0, 0])},
        "S2": {"G": _sat_day([0, 0, 0, 1, 1, 1, 1])},
    }
    plain = ps.simulate_constellation([day], step_s=60, info_bitrate_bps=1e6, generated_bps=1e4)
    setup = ps.simulate_constellation([day], step_s=60, info_bitrate_bps=1e6, generated_bps=1e4, setup_s=120)
    assert plain.satellites[1]["served_min_per_day"] == 4.0
    assert setup.satellites[1]["served_min_per_day"] == 2.0  # samples 5-6
    assert setup.satellites[1]["conflict_min_per_day"] == 2.0  # waiting for the antenna
    assert setup.satellites[0]["served_min_per_day"] == 2.0


def test_statistics_report_pdf(tmp_path):
    import statistics_report

    d1 = _day(1, [0, 1, 1, 0, 1, 0], [0, 4.0, 1.0, 0, 5.0, 0], station=["", "G", "G", "", "G", ""])
    result = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, min_margin_db=2.0, generated_bps=1e4,
                          latency_requirement=(45.0, 90.0))
    path = tmp_path / "report.pdf"
    statistics_report.build_statistics_report_pdf(
        str(path), result,
        statistics_report.StatisticsReportContext(
            title="TEST – Statistics", period="TEST over G, 1 day", settings=[("Satellite", "TEST")],
            summary_lines=["Passes: 2", "Usable time (margin ≥ 2 dB): 20 min/day"], project="Demo",
            app_version="9.9.9",
        ),
    )
    data = path.read_bytes()
    assert data.startswith(b"%PDF") and len(data) > 2000


def test_drop_short_contacts_across_midnight():
    start1 = datetime.datetime(2026, 10, 1, tzinfo=UTC)
    start2 = datetime.datetime(2026, 10, 2, tzinfo=UTC)

    def samples(start, visible):
        visible = np.array(visible, dtype=bool)
        return ps.DaySamples(start, np.where(visible, 20.0, -5.0), visible, np.where(visible, 5.0, np.nan))

    # Contact 1: one sample (dropped); contact 2: two samples before and one after midnight (kept).
    days = [
        {"S": {"G": samples(start1, [0, 1, 0, 0, 1, 1])}},
        {"S": {"G": samples(start2, [1, 0, 0, 0, 0, 0])}},
    ]
    kept = ps.drop_short_contacts(days, step_s=60, min_contact_s=120)
    assert list(kept[0]["S"]["G"].visible) == [False, False, False, False, True, True]
    assert list(kept[1]["S"]["G"].visible) == [True, False, False, False, False, False]
    assert days[0]["S"]["G"].visible[1]  # the input is not changed
    assert ps.drop_short_contacts(days, 60, 0) is days


def test_latency_at_pass_end():
    # Contact in the last two samples: with the dump counted at its end, every bit waits one sample more
    # when it came down in the first sample of the pass.
    d1 = _day(1, [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 5, 5])
    per_bit = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4)
    at_end = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4, latency_at_pass_end=True)
    assert at_end.summary["max_latency_min"] == pytest.approx(60.0)  # all on ground at the end of the pass
    assert at_end.summary["mean_latency_min"] > per_bit.summary["mean_latency_min"]
    day = {"S1": {"G": _sat_day([0, 0, 0, 0, 1, 1])}}
    c_bit = ps.simulate_constellation([day], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4)
    c_end = ps.simulate_constellation([day], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4,
                                      latency_at_pass_end=True)
    assert c_end.summary["max_latency_min"] == pytest.approx(60.0)
    assert c_end.summary["mean_latency_min"] > c_bit.summary["mean_latency_min"]


def test_memory_capacity_and_loss():
    # 0.75 MB generated per 10-minute sample, contact only in the last sample: up to 4.5 MB on board.
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    roomy = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4, memory_bits=8 * 9e6)
    assert roomy.summary["max_memory_use_pct"] == pytest.approx(3.75 / 9 * 100)  # 3.75 MB before the dump
    assert roomy.summary["lost_mb"] == 0.0
    assert roomy.summary["memory_autonomy_h"] == pytest.approx(9e6 * 8 / 1e4 / 3600)
    small = ps.aggregate([d1], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4, memory_bits=8 * 2e6)
    assert small.summary["lost_mb"] == pytest.approx(3.75 - 2.0 + 0.0, abs=0.76)
    assert small.summary["max_memory_use_pct"] == pytest.approx(100.0)
    day = {"S": {"G": _sat_day([0, 0, 0, 0, 0, 1])}}
    c = ps.simulate_constellation([day], step_s=600, info_bitrate_bps=1e6, generated_bps=1e4, memory_bits=8 * 2e6)
    assert c.summary["lost_mb"] > 0 and c.satellites[0]["lost_mb"] > 0


def test_products_on_ground():
    # Data of 0-30 min down at 35 min, of 30-60 min at 70 min, of 60-75 min at 80 min; products of 60 min.
    raw = [(1800.0, 0.0, 1800.0, 2100.0), (1800.0, 1800.0, 3600.0, 4200.0), (900.0, 3600.0, 4500.0, 4800.0)]
    data = ps._chunks(raw, 3600.0, "data")
    # The first product is complete at 70 min; the second (60-120 min) never is and is not counted.
    assert sum(c[0] for c in data) == pytest.approx(3600.0)
    assert min(c[1] for c in data) == pytest.approx(600.0)  # newest datum: 70 - 60 min
    assert max(c[2] for c in data) == pytest.approx(4200.0)  # oldest datum: 70 - 0 min
    assert [c[1] for c in ps._chunks(raw, 3600.0, "start")] == pytest.approx([4200.0, 4200.0])
    assert [c[1] for c in ps._chunks(raw, 3600.0, "end")] == pytest.approx([600.0, 600.0])
    # Without products every chunk is on ground when it comes down.
    assert ps._chunks(raw)[0] == (1800.0, 300.0, 2100.0)
    with pytest.raises(ValueError):
        ps._chunks(raw, 3600.0, "middle")


def test_products_in_the_statistics():
    # One dump at the end of the hour: two 30-min products, both complete at 60 min.
    d1 = _day(1, [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 5])
    kwargs = dict(step_s=600, info_bitrate_bps=1e6, generated_bps=1e4)
    per_bit = ps.aggregate([d1], **kwargs)
    by_data = ps.aggregate([d1], product_s=1800, **kwargs)
    by_start = ps.aggregate([d1], product_s=1800, product_reference="start", **kwargs)
    by_end = ps.aggregate([d1], product_s=1800, product_reference="end", **kwargs)
    assert by_data.summary["mean_latency_min"] == pytest.approx(per_bit.summary["mean_latency_min"])
    assert by_start.summary["mean_latency_min"] == pytest.approx(45.0)  # 60 and 30 min
    assert by_end.summary["mean_latency_min"] == pytest.approx(15.0)  # 30 and 0 min
    day = {"S1": {"G": _sat_day([0, 0, 0, 0, 0, 1])}}
    c_start = ps.simulate_constellation([day], product_s=1800, product_reference="start", **kwargs)
    assert c_start.summary["mean_latency_min"] == pytest.approx(45.0)
