"""Contact and data volume statistics over several days and ground stations.

The pass analysis looks at one day and one station; this module adds up
its samples over a period and a network of stations: the contacts of each
day, the minutes in which the link closes (margin at least the required
one), the gaps between contacts and the data downloaded.

Without on-board data the volume is the downlink capacity: the information
bit rate over the usable minutes. With the data generated on board (the
instruments, per orbit, and the housekeeping telemetry, HKTM) the on-board
memory is simulated: it fills at the generation rate and empties at the
downlink rate during the usable minutes, never below zero. This gives the
volume actually downloaded, the memory occupation and the data latency
(time from generation to arrival on ground, first in first out).

It only aggregates: the samples (elevation, whether the satellite is in
contact, link margin) are computed by the caller, one day and one station
at a time, with the same models as the pass analysis.
"""
from __future__ import annotations

import datetime
from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional

import numpy as np


@dataclass
class DaySamples:
    """Samples of one day, every ``step_s`` seconds from ``start``."""

    start: datetime.datetime
    elevation_deg: np.ndarray
    visible: np.ndarray
    # Link margin (dB) over the required one; None for a contact analysis only.
    margin_db: Optional[np.ndarray] = None
    # Station of each sample in contact (the best one of a network), else "".
    station: Optional[np.ndarray] = None


@dataclass
class PassRecord:
    start: datetime.datetime
    end: datetime.datetime
    max_elevation_deg: float
    contact_min: float
    stations: str = ""
    usable_min: Optional[float] = None
    capacity_mb: Optional[float] = None
    downloaded_mb: Optional[float] = None
    min_margin_db: Optional[float] = None


@dataclass
class DayRecord:
    date: datetime.date
    contacts: int = 0
    contact_min: float = 0.0
    usable_min: Optional[float] = None
    capacity_mb: Optional[float] = None
    downloaded_mb: Optional[float] = None
    max_gap_min: Optional[float] = None
    on_board_mb: Optional[float] = None


@dataclass
class StatisticsResult:
    passes: list[PassRecord] = field(default_factory=list)
    days: list[DayRecord] = field(default_factory=list)
    summary: dict = field(default_factory=dict)


def combine_stations(samples: dict[str, DaySamples]) -> DaySamples:
    """One day of a station network: in contact when any station is.

    At each sample the station with the best margin (or, for a contact
    analysis, the highest elevation) is the one used.
    """
    names = list(samples)
    first = samples[names[0]]
    visible = np.vstack([np.asarray(samples[n].visible, dtype=bool) for n in names])
    elevation = np.vstack([samples[n].elevation_deg for n in names])
    with_budget = first.margin_db is not None
    if with_budget:
        score = np.vstack([np.where(samples[n].visible, samples[n].margin_db, -np.inf) for n in names])
    else:
        score = np.where(visible, elevation, -np.inf)
    best = np.argmax(score, axis=0)
    columns = np.arange(visible.shape[1])
    any_visible = visible.any(axis=0)
    margin = None
    if with_budget:
        margin_all = np.vstack([samples[n].margin_db for n in names])
        margin = np.where(any_visible, margin_all[best, columns], np.nan)
    station = np.where(any_visible, np.array(names, dtype=object)[best], "")
    return DaySamples(
        first.start,
        np.where(any_visible, elevation[best, columns], elevation.max(axis=0)),
        any_visible,
        margin,
        station,
    )


def _groups(flags: np.ndarray) -> list[tuple[int, int]]:
    """(first, last) indices of each run of True in ``flags``."""
    runs = []
    start = None
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(flags) - 1))
    return runs


def _ordered_unique(values) -> list:
    seen = []
    for value in values:
        if value and value not in seen:
            seen.append(value)
    return seen


def aggregate(
    days: Iterable[DaySamples],
    step_s: float,
    info_bitrate_bps: Optional[float] = None,
    min_margin_db: float = 0.0,
    generated_bps: Optional[float] = None,
) -> StatisticsResult:
    """Passes, daily figures and a summary from the samples of consecutive days.

    A pass is a run of samples in contact; one crossing midnight counts for
    the day it starts. Its usable minutes are those with ``margin_db >=
    min_margin_db``; the capacity is ``info_bitrate_bps`` over them (MB of
    10^6 bytes). With ``generated_bps`` (average rate of the data produced
    on board) the memory is simulated, starting empty: the downloaded volume
    is at most what is on board. The gap after a pass lasts until the next
    pass starts and counts for the day the next pass starts.
    """
    days = list(days)
    result = StatisticsResult()
    if not days:
        return result
    step_min = step_s / 60.0
    with_budget = days[0].margin_db is not None and info_bitrate_bps is not None
    with_memory = with_budget and generated_bps is not None and generated_bps > 0

    elevation = np.concatenate([d.elevation_deg for d in days])
    visible = np.concatenate([np.asarray(d.visible, dtype=bool) for d in days])
    n = len(visible)
    station = (
        np.concatenate([d.station for d in days]) if days[0].station is not None else np.full(n, "", dtype=object)
    )
    times = [d.start + datetime.timedelta(seconds=step_s * i) for d in days for i in range(len(d.elevation_deg))]
    day_index = np.concatenate([np.full(len(d.elevation_deg), k) for k, d in enumerate(days)])

    capacity_bits = np.zeros(n)
    downloaded_bits = np.zeros(n)
    margin = None
    if with_budget:
        margin = np.concatenate([d.margin_db for d in days])
        usable = visible & (np.nan_to_num(margin, nan=-np.inf) >= min_margin_db)
        capacity_bits = np.where(usable, info_bitrate_bps * step_s, 0.0)
        downloaded_bits = capacity_bits.copy()

    memory_end_bits = None
    latency = {}
    if with_memory:
        memory = 0.0
        memory_end_bits = np.zeros(n)
        downloaded_total = 0.0
        latencies_weighted = 0.0
        max_latency = 0.0
        generated_per_step = generated_bps * step_s
        for i in range(n):
            memory += generated_per_step
            take = min(memory, capacity_bits[i])
            if take > 0:
                # First in, first out: the bits downloaded now were generated
                # from (downloaded_total) to (downloaded_total + take) / rate.
                t_now = (i + 1) * step_s
                oldest = downloaded_total / generated_bps
                newest = (downloaded_total + take) / generated_bps
                max_latency = max(max_latency, t_now - oldest)
                latencies_weighted += take * (t_now - (oldest + newest) / 2)
                downloaded_total += take
            memory -= take
            downloaded_bits[i] = take
            memory_end_bits[i] = memory
        latency = {
            "max_latency_min": max_latency / 60.0 if downloaded_total else None,
            "mean_latency_min": latencies_weighted / downloaded_total / 60.0 if downloaded_total else None,
        }

    by_date = [DayRecord(d.start.date()) for d in days]
    for k, record in enumerate(by_date):
        in_day = day_index == k
        if with_budget:
            record.usable_min = float(np.count_nonzero(capacity_bits[in_day])) * step_min
            record.capacity_mb = float(capacity_bits[in_day].sum()) / 8e6
            record.downloaded_mb = float(downloaded_bits[in_day].sum()) / 8e6
        if memory_end_bits is not None:
            record.on_board_mb = float(memory_end_bits[in_day][-1]) / 8e6

    previous_end = None
    gaps = []
    for first, last in _groups(visible):
        record = PassRecord(
            start=times[first],
            end=times[last] + datetime.timedelta(seconds=step_s),
            max_elevation_deg=float(elevation[first : last + 1].max()),
            contact_min=(last - first + 1) * step_min,
            stations=" + ".join(_ordered_unique(station[first : last + 1])),
        )
        day = by_date[int(day_index[first])]
        day.contacts += 1
        day.contact_min += record.contact_min
        if with_budget:
            record.usable_min = float(np.count_nonzero(capacity_bits[first : last + 1])) * step_min
            record.capacity_mb = float(capacity_bits[first : last + 1].sum()) / 8e6
            record.downloaded_mb = float(downloaded_bits[first : last + 1].sum()) / 8e6
            record.min_margin_db = float(np.nanmin(margin[first : last + 1]))
        if previous_end is not None:
            gap = (record.start - previous_end).total_seconds() / 60.0
            gaps.append(gap)
            day.max_gap_min = gap if day.max_gap_min is None else max(day.max_gap_min, gap)
        previous_end = record.end
        result.passes.append(record)
    result.days = by_date

    n_days = len(by_date)
    summary = {
        "days": n_days,
        "passes": len(result.passes),
        "passes_per_day": len(result.passes) / n_days,
        "contact_min_per_day": sum(d.contact_min for d in by_date) / n_days,
        "max_gap_min": max(gaps) if gaps else None,
        "mean_gap_min": float(np.mean(gaps)) if gaps else None,
    }
    if with_budget:
        capacity = [d.capacity_mb for d in by_date]
        downloaded = [d.downloaded_mb for d in by_date]
        summary.update(
            usable_min_per_day=sum(d.usable_min for d in by_date) / n_days,
            capacity_mb_per_day_mean=float(np.mean(capacity)),
            capacity_mb_per_day_min=min(capacity),
            capacity_mb_per_day_max=max(capacity),
            downloaded_mb_total=sum(downloaded),
            downloaded_mb_per_day_mean=float(np.mean(downloaded)),
            info_bitrate_bps=info_bitrate_bps,
            min_margin_db=min_margin_db,
        )
    if with_memory:
        summary.update(
            generated_bps=generated_bps,
            generated_mb_per_day=generated_bps * 86400 / 8e6,
            max_on_board_mb=float(memory_end_bits.max()) / 8e6,
            final_on_board_mb=float(memory_end_bits[-1]) / 8e6,
            **latency,
        )
    result.summary = summary
    return result


def compute(
    day_samples: Callable[[datetime.datetime], DaySamples],
    start_date: datetime.date,
    n_days: int,
    step_s: float,
    progress: Optional[Callable[[int, int], None]] = None,
    **aggregate_kwargs,
) -> StatisticsResult:
    """Run ``day_samples`` for each day from ``start_date`` and aggregate the results."""
    if n_days < 1:
        raise ValueError("The number of days must be at least 1.")
    samples = []
    for k in range(n_days):
        day = datetime.datetime.combine(
            start_date + datetime.timedelta(days=k), datetime.time(0, 0), tzinfo=datetime.timezone.utc
        )
        samples.append(day_samples(day))
        if progress is not None:
            progress(k + 1, n_days)
    return aggregate(samples, step_s, **aggregate_kwargs)
