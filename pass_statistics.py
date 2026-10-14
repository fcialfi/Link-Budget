"""Contact and data volume statistics over several days.

The pass analysis looks at one day; this module adds up its samples over a
period: the contacts of each day, the minutes in which the link closes
(margin at least the required one), the data volume downloaded in those
minutes and the gaps between contacts.

It only aggregates: the samples (elevation, whether the satellite is in
contact, link margin) are computed by the caller, one day at a time, with
the same models as the pass analysis.
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


@dataclass
class PassRecord:
    start: datetime.datetime
    end: datetime.datetime
    max_elevation_deg: float
    contact_min: float
    usable_min: Optional[float] = None
    volume_mb: Optional[float] = None
    min_margin_db: Optional[float] = None


@dataclass
class DayRecord:
    date: datetime.date
    contacts: int = 0
    contact_min: float = 0.0
    usable_min: Optional[float] = None
    volume_mb: Optional[float] = None
    max_gap_min: Optional[float] = None
    backlog_mb: Optional[float] = None


@dataclass
class StatisticsResult:
    passes: list[PassRecord] = field(default_factory=list)
    days: list[DayRecord] = field(default_factory=list)
    summary: dict = field(default_factory=dict)


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


def aggregate(
    days: Iterable[DaySamples],
    step_s: float,
    info_bitrate_bps: Optional[float] = None,
    min_margin_db: float = 0.0,
    onboard_mb_per_day: Optional[float] = None,
) -> StatisticsResult:
    """Passes, daily figures and a summary from the samples of consecutive days.

    A pass is a run of samples in contact; one crossing midnight counts for
    the day it starts. Its usable minutes are those with ``margin_db >=
    min_margin_db``; the volume is ``info_bitrate_bps`` over them (MB of
    10^6 bytes). The gap after a pass lasts until the next pass starts and
    counts for the day the next pass starts. With ``onboard_mb_per_day`` the
    backlog at the end of each day is what is generated on board and not yet
    downloaded (never below 0).
    """
    days = list(days)
    result = StatisticsResult()
    if not days:
        return result
    step_min = step_s / 60.0
    with_budget = days[0].margin_db is not None and info_bitrate_bps is not None

    # One timeline: the samples of all days back to back.
    elevation = np.concatenate([d.elevation_deg for d in days])
    visible = np.concatenate([np.asarray(d.visible, dtype=bool) for d in days])
    margin = np.concatenate([d.margin_db for d in days]) if with_budget else None
    times = [d.start + datetime.timedelta(seconds=step_s * i) for d in days for i in range(len(d.elevation_deg))]

    by_date = {d.start.date(): DayRecord(d.start.date()) for d in days}
    for record in by_date.values():
        if with_budget:
            record.usable_min = 0.0
            record.volume_mb = 0.0

    previous_end = None
    gaps = []
    for first, last in _groups(visible):
        record = PassRecord(
            start=times[first],
            end=times[last] + datetime.timedelta(seconds=step_s),
            max_elevation_deg=float(elevation[first : last + 1].max()),
            contact_min=(last - first + 1) * step_min,
        )
        day = by_date[record.start.date()]
        day.contacts += 1
        day.contact_min += record.contact_min
        if with_budget:
            pass_margin = margin[first : last + 1]
            usable = int(np.count_nonzero(pass_margin >= min_margin_db))
            record.usable_min = usable * step_min
            record.volume_mb = usable * step_s * info_bitrate_bps / 8e6
            record.min_margin_db = float(np.nanmin(pass_margin))
            day.usable_min += record.usable_min
            day.volume_mb += record.volume_mb
        if previous_end is not None:
            gap = (record.start - previous_end).total_seconds() / 60.0
            gaps.append(gap)
            day.max_gap_min = gap if day.max_gap_min is None else max(day.max_gap_min, gap)
        previous_end = record.end
        result.passes.append(record)

    result.days = [by_date[d.start.date()] for d in days]
    if with_budget and onboard_mb_per_day is not None:
        backlog = 0.0
        for day in result.days:
            backlog = max(0.0, backlog + onboard_mb_per_day - day.volume_mb)
            day.backlog_mb = backlog

    n_days = len(result.days)
    summary = {
        "days": n_days,
        "passes": len(result.passes),
        "passes_per_day": len(result.passes) / n_days,
        "contact_min_per_day": sum(d.contact_min for d in result.days) / n_days,
        "max_gap_min": max(gaps) if gaps else None,
        "mean_gap_min": float(np.mean(gaps)) if gaps else None,
    }
    if with_budget:
        volumes = [d.volume_mb for d in result.days]
        summary.update(
            usable_min_per_day=sum(d.usable_min for d in result.days) / n_days,
            volume_mb_total=sum(volumes),
            volume_mb_per_day_mean=float(np.mean(volumes)),
            volume_mb_per_day_min=min(volumes),
            volume_mb_per_day_max=max(volumes),
            info_bitrate_bps=info_bitrate_bps,
            min_margin_db=min_margin_db,
        )
        if onboard_mb_per_day is not None:
            summary["onboard_mb_per_day"] = onboard_mb_per_day
            summary["final_backlog_mb"] = result.days[-1].backlog_mb
            summary["max_backlog_mb"] = max(d.backlog_mb for d in result.days)
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
