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
    # Downloaded chunks: (bits, smallest latency s, largest latency s); the bits
    # of a chunk were generated uniformly over that latency range.
    latency_chunks: list[tuple[float, float, float]] = field(default_factory=list)

    def latency_fraction(self, latency_min: float) -> Optional[float]:
        """Share (0-1) of the downloaded data that reached the ground within ``latency_min``."""
        total = sum(bits for bits, _lo, _hi in self.latency_chunks)
        if not total:
            return None
        limit = latency_min * 60.0
        within = 0.0
        for bits, lo, hi in self.latency_chunks:
            if hi <= limit:
                within += bits
            elif lo < limit:
                within += bits * (limit - lo) / (hi - lo)
        return within / total

    def latency_percentile(self, percent: float) -> Optional[float]:
        """Latency (min) within which ``percent`` % of the downloaded data reached the ground."""
        if not self.latency_chunks:
            return None
        target = percent / 100.0
        low, high = 0.0, max(hi for _bits, _lo, hi in self.latency_chunks) / 60.0
        for _ in range(60):
            mid = (low + high) / 2
            if self.latency_fraction(mid) >= target:
                high = mid
            else:
                low = mid
        return high


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
    latency_requirement: Optional[tuple[float, float]] = None,
) -> StatisticsResult:
    """Passes, daily figures and a summary from the samples of consecutive days.

    A pass is a run of samples in contact; one crossing midnight counts for
    the day it starts. Its usable minutes are those with ``margin_db >=
    min_margin_db``; the capacity is ``info_bitrate_bps`` over them (MB of
    10^6 bytes). With ``generated_bps`` (average rate of the data produced
    on board) the memory is simulated, starting empty: the downloaded volume
    is at most what is on board. The gap after a pass lasts until the next
    pass starts and counts for the day the next pass starts.

    ``latency_requirement`` = (minutes, percent), e.g. (60, 90) for "90 % of
    the data on ground within 60 min", is checked on the downloaded data.
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
                result.latency_chunks.append((take, t_now - newest, t_now - oldest))
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
        if latency_requirement is not None and result.latency_chunks:
            minutes, percent = latency_requirement
            within = float(result.latency_fraction(minutes)) * 100.0
            summary.update(
                latency_requirement_min=minutes,
                latency_requirement_pct=percent,
                latency_within_pct=within,
                latency_at_required_pct_min=result.latency_percentile(percent),
                latency_requirement_met=bool(within >= percent),
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


# ---------------------------------------------------------------------------
# Constellation: several satellites sharing the ground stations
# ---------------------------------------------------------------------------


@dataclass
class ConstellationResult(StatisticsResult):
    # One dict per satellite and per station (see simulate_constellation).
    satellites: list[dict] = field(default_factory=list)
    stations: list[dict] = field(default_factory=list)


def simulate_constellation(
    days: list[dict[str, dict[str, DaySamples]]],
    step_s: float,
    info_bitrate_bps: float,
    min_margin_db: float = 0.0,
    generated_bps: float = 0.0,
    latency_requirement: Optional[tuple[float, float]] = None,
) -> ConstellationResult:
    """Satellites sharing ground stations that track one satellite at a time.

    ``days[k][satellite][station]`` are the samples of day k. Every satellite
    generates ``generated_bps`` into its own memory (empty at the start). At
    each sample each station serves at most one of the satellites for which
    the link closes (margin >= ``min_margin_db``), and each satellite uses at
    most one station: the satellites with more data on board are served
    first, each through its best-margin station still free. A satellite that
    could have downloaded but found its stations busy loses that sample
    (conflict). Latency is first in first out per satellite; the system
    figures put the data of all satellites together.
    """
    result = ConstellationResult()
    if not days:
        return result
    sat_names = list(days[0])
    station_names = list(days[0][sat_names[0]])
    n_days = len(days)
    per_day = [len(days[k][sat_names[0]][station_names[0]].visible) for k in range(n_days)]
    n = sum(per_day)
    starts = [days[k][sat_names[0]][station_names[0]].start for k in range(n_days)]
    times = [starts[k] + datetime.timedelta(seconds=step_s * i) for k in range(n_days) for i in range(per_day[k])]
    day_index = np.concatenate([np.full(per_day[k], k) for k in range(n_days)])

    def joined(sat, st, attr):
        return np.concatenate([getattr(days[k][sat][st], attr) for k in range(n_days)])

    # usable[s, g, i], margin[s, g, i], elevation[s, g, i]
    margin = np.stack([np.stack([joined(s, g, "margin_db") for g in station_names]) for s in sat_names])
    visible = np.stack([np.stack([joined(s, g, "visible") for g in station_names]) for s in sat_names]).astype(bool)
    elevation = np.stack([np.stack([joined(s, g, "elevation_deg") for g in station_names]) for s in sat_names])
    usable = visible & (np.nan_to_num(margin, nan=-np.inf) >= min_margin_db)

    n_sat, n_st = len(sat_names), len(station_names)
    step_bits = info_bitrate_bps * step_s
    gen_step = generated_bps * step_s
    memory = np.zeros(n_sat)
    downloaded_total = np.zeros(n_sat)
    served_station = np.full((n_sat, n), -1)
    downloaded = np.zeros((n_sat, n))
    on_board = np.zeros((n_sat, n))
    conflict = np.zeros((n_sat, n), dtype=bool)
    chunks: list[list] = [[] for _ in range(n_sat)]
    for i in range(n):
        memory += gen_step
        # A station keeps tracking its satellite while the link closes (a real
        # antenna follows one pass); the free stations then go to the other
        # satellites by a maximum matching (augmenting paths), the satellites
        # with more data on board first.
        owner: dict[int, int] = {}
        if i > 0:
            for s in range(n_sat):
                g = served_station[s, i - 1]
                if g >= 0 and usable[s, g, i]:
                    owner[int(g)] = s
        tracked = set(owner.values())
        options = {
            s: sorted(
                (g for g in np.flatnonzero(usable[s, :, i]) if int(g) not in owner),
                key=lambda g, s=s: -margin[s, g, i],
            )
            for s in range(n_sat)
        }
        locked = set(owner)

        def _assign(s, seen):
            for g in options[s]:
                g = int(g)
                if g in seen or g in locked:
                    continue
                seen.add(g)
                if g not in owner or _assign(owner[g], seen):
                    owner[g] = s
                    return True
            return False

        for s in np.argsort(-memory, kind="stable"):
            s = int(s)
            if s in tracked:
                continue
            if usable[s, :, i].any() and not _assign(s, set()):
                conflict[s, i] = True
        for g, s in owner.items():
            served_station[s, i] = g
        for s in range(n_sat):
            if served_station[s, i] < 0:
                continue
            take = min(memory[s], step_bits)
            if take > 0 and generated_bps > 0:
                t_now = (i + 1) * step_s
                oldest = downloaded_total[s] / generated_bps
                newest = (downloaded_total[s] + take) / generated_bps
                chunks[s].append((take, t_now - newest, t_now - oldest))
            elif generated_bps <= 0:
                take = step_bits  # capacity only
            memory[s] -= take if generated_bps > 0 else 0.0
            downloaded_total[s] += take
            downloaded[s, i] = take
        on_board[:, i] = memory

    step_min = step_s / 60.0
    with_memory = generated_bps > 0
    for s, name in enumerate(sat_names):
        sat_chunks = chunks[s]
        sat_result = StatisticsResult(latency_chunks=sat_chunks)
        bits = sum(c[0] for c in sat_chunks)
        record = {
            "satellite": name,
            "served_min_per_day": float(np.count_nonzero(served_station[s] >= 0)) * step_min / n_days,
            "conflict_min_per_day": float(np.count_nonzero(conflict[s])) * step_min / n_days,
            "downloaded_mb_per_day": float(downloaded[s].sum()) / 8e6 / n_days,
            "max_on_board_mb": float(on_board[s].max()) / 8e6 if with_memory else None,
            "final_on_board_mb": float(on_board[s, -1]) / 8e6 if with_memory else None,
            "mean_latency_min": (
                sum(c[0] * (c[1] + c[2]) / 2 for c in sat_chunks) / bits / 60.0 if bits else None
            ),
            "max_latency_min": max(c[2] for c in sat_chunks) / 60.0 if sat_chunks else None,
        }
        if latency_requirement is not None and sat_chunks:
            record["latency_within_pct"] = float(sat_result.latency_fraction(latency_requirement[0])) * 100.0
        result.satellites.append(record)
        # Service windows of this satellite: runs of served samples.
        for first, last in _groups(served_station[s] >= 0):
            used = [station_names[g] for g in served_station[s, first : last + 1]]
            best_el = [elevation[s, served_station[s, i], i] for i in range(first, last + 1)]
            served_margin = [margin[s, served_station[s, i], i] for i in range(first, last + 1)]
            result.passes.append(
                PassRecord(
                    start=times[first],
                    end=times[last] + datetime.timedelta(seconds=step_s),
                    max_elevation_deg=float(max(best_el)),
                    contact_min=(last - first + 1) * step_min,
                    stations=f"{name}: " + " + ".join(_ordered_unique(used)),
                    usable_min=(last - first + 1) * step_min,
                    capacity_mb=(last - first + 1) * step_bits / 8e6,
                    downloaded_mb=float(downloaded[s, first : last + 1].sum()) / 8e6,
                    min_margin_db=float(np.nanmin(served_margin)),
                )
            )
    result.passes.sort(key=lambda p: p.start)

    for g, name in enumerate(station_names):
        served_here = served_station == g
        result.stations.append({
            "station": name,
            "served_min_per_day": float(np.count_nonzero(served_here)) * step_min / n_days,
            "downloaded_mb_per_day": float(downloaded[served_here].sum()) / 8e6 / n_days,
            "busy_pct": float(np.count_nonzero(served_here.any(axis=0))) / n * 100.0,
            "contention_min_per_day": float(
                np.count_nonzero(usable[:, g, :].sum(axis=0) > 1)
            ) * step_min / n_days,
        })

    for k in range(n_days):
        in_day = day_index == k
        record = DayRecord(starts[k].date())
        record.contacts = sum(1 for p in result.passes if p.start.date() == record.date)
        record.contact_min = float(np.count_nonzero(served_station[:, in_day] >= 0)) * step_min
        record.usable_min = float(np.count_nonzero(usable[:, :, in_day].any(axis=1))) * step_min
        record.capacity_mb = float(np.count_nonzero(usable[:, :, in_day].any(axis=1))) * step_bits / 8e6
        record.downloaded_mb = float(downloaded[:, in_day].sum()) / 8e6
        if with_memory:
            record.on_board_mb = float(on_board[:, in_day][:, -1].sum()) / 8e6
        result.days.append(record)

    all_chunks = [c for sat_chunks in chunks for c in sat_chunks]
    result.latency_chunks = all_chunks
    total_bits = sum(c[0] for c in all_chunks)
    summary = {
        "satellites": n_sat,
        "stations": n_st,
        "days": n_days,
        "served_min_per_day": float(np.count_nonzero(served_station >= 0)) * step_min / n_days,
        "conflict_min_per_day": float(np.count_nonzero(conflict)) * step_min / n_days,
        "downloaded_mb_per_day": float(downloaded.sum()) / 8e6 / n_days,
        "info_bitrate_bps": info_bitrate_bps,
        "min_margin_db": min_margin_db,
    }
    if with_memory:
        summary.update(
            generated_bps=generated_bps,
            generated_mb_per_day=generated_bps * 86400 / 8e6 * n_sat,
            max_on_board_mb=float(on_board.sum(axis=0).max()) / 8e6,
            final_on_board_mb=float(on_board[:, -1].sum()) / 8e6,
            mean_latency_min=sum(c[0] * (c[1] + c[2]) / 2 for c in all_chunks) / total_bits / 60.0 if total_bits else None,
            max_latency_min=max(c[2] for c in all_chunks) / 60.0 if all_chunks else None,
        )
        if latency_requirement is not None and all_chunks:
            minutes, percent = latency_requirement
            within = float(result.latency_fraction(minutes)) * 100.0
            summary.update(
                latency_requirement_min=minutes,
                latency_requirement_pct=percent,
                latency_within_pct=within,
                latency_at_required_pct_min=result.latency_percentile(percent),
                latency_requirement_met=bool(within >= percent),
            )
    result.summary = summary
    return result
