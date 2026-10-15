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
    # Data on board (MB) at the end of every sample, with the sample times; for a
    # constellation, the total of all satellites. None without the memory model.
    on_board_series: Optional[tuple[list, np.ndarray]] = None

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
            elif lo < limit:  # here hi > limit > lo
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


def _memory_figures(memory_bits, max_bits_one_satellite, lost_bits, generated_bps, period_s) -> dict:
    """Use of the on-board memory of one satellite: share, autonomy without contacts, data lost."""
    if not memory_bits:
        return {}
    return {
        "memory_mb": memory_bits / 8e6,
        "max_memory_use_pct": 100.0 * max_bits_one_satellite / memory_bits,
        "memory_autonomy_h": memory_bits / generated_bps / 3600.0,
        "lost_mb": lost_bits / 8e6,
        "lost_pct": 100.0 * lost_bits / (generated_bps * period_s) if period_s else 0.0,
    }


PRODUCT_REFERENCES = ("data", "start", "end")


def _chunks(raw, product_s: Optional[float] = None, reference: str = "data") -> list[tuple[float, float, float]]:
    """Latency chunks (bits, smallest latency s, largest latency s) of one satellite.

    ``raw`` lists what came down, first in first out: (bits, generated from s,
    generated to s, on ground at s), times from the start of the period. Without
    ``product_s`` each chunk is on ground when it comes down. With it the ground
    builds products of ``product_s`` seconds of on-board data, fixed blocks from
    the start of the period (00:00 UTC): a product is complete when the last of
    its data is on ground, and only complete products are counted. Its data then
    have the latency from ``reference``: "data", the sensing time of each datum;
    "start" or "end", the start or end of the product's sensing.
    """
    if not product_s:
        return [(bits, t_arr - g1, t_arr - g0) for bits, g0, g1, t_arr in raw]
    if reference not in PRODUCT_REFERENCES:
        raise ValueError(f"Unknown product latency reference '{reference}' (use data, start or end).")
    pieces = []  # (product, bits, generated from, generated to)
    arrival: dict[int, float] = {}
    covered = 0.0  # generation time downloaded so far (first in first out)
    for bits, g0, g1, t_arr in raw:
        if bits <= 0 or g1 <= g0:
            continue
        covered = max(covered, g1)
        k = int(g0 // product_s)
        while True:
            lo, hi = max(g0, k * product_s), min(g1, (k + 1) * product_s)
            if hi > lo:
                pieces.append((k, bits * (hi - lo) / (g1 - g0), lo, hi))
                arrival[k] = max(arrival.get(k, 0.0), t_arr)
            if g1 <= (k + 1) * product_s:
                break
            k += 1
    result = []
    for k, bits, lo, hi in pieces:
        if covered < (k + 1) * product_s - 1e-6:
            continue  # product not complete at the end of the period
        t_arr = arrival[k]
        if reference == "data":
            result.append((bits, t_arr - hi, t_arr - lo))
        else:
            latency = t_arr - (k * product_s if reference == "start" else (k + 1) * product_s)
            result.append((bits, latency, latency))
    return result


def _latency_figures(chunks) -> dict:
    total = sum(c[0] for c in chunks)
    return {
        "max_latency_min": max(c[2] for c in chunks) / 60.0 if chunks else None,
        "mean_latency_min": sum(c[0] * (c[1] + c[2]) / 2 for c in chunks) / total / 60.0 if total else None,
    }


def drop_short_contacts(
    days: list[dict[str, dict[str, DaySamples]]],
    step_s: float,
    min_contact_s: float,
    min_margin_db: float = 0.0,
) -> list[dict[str, dict[str, DaySamples]]]:
    """The same days without the contacts shorter than ``min_contact_s``.

    A contact is a run of usable samples (in view and, with a link budget,
    margin >= ``min_margin_db``) of one satellite over one station, followed
    across midnight. The samples of a shorter one are set out of view: a pass
    too short to be worth scheduling (acquisition, lock, dump start).
    """
    if not days or min_contact_s <= 0:
        return days
    min_samples = int(np.ceil(min_contact_s / step_s))
    result = [{sat: dict(per_station) for sat, per_station in day.items()} for day in days]
    for sat in days[0]:
        for station in days[0][sat]:
            parts = [day[sat][station] for day in days]
            usable = np.concatenate([
                np.asarray(p.visible, dtype=bool)
                & (True if p.margin_db is None else np.nan_to_num(p.margin_db, nan=-np.inf) >= min_margin_db)
                for p in parts
            ])
            keep = np.ones(len(usable), dtype=bool)
            for first, last in _groups(usable):
                if last - first + 1 < min_samples:
                    keep[first : last + 1] = False
            offset = 0
            for k, p in enumerate(parts):
                n = len(p.visible)
                part_keep = keep[offset : offset + n]
                offset += n
                if part_keep.all():
                    continue
                margin = None if p.margin_db is None else np.where(part_keep, p.margin_db, np.nan)
                result[k][sat][station] = DaySamples(
                    p.start, p.elevation_deg, np.asarray(p.visible, dtype=bool) & part_keep, margin, p.station
                )
    return result


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
    latency_at_pass_end: bool = False,
    memory_bits: Optional[float] = None,
    product_s: Optional[float] = None,
    product_reference: str = "data",
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
    ``product_s`` groups the data on ground into products of that many seconds
    of on-board data (see :func:`_chunks`).
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
        raw = []
        generated_per_step = generated_bps * step_s
        # With latency_at_pass_end the data count as on ground at the end of
        # their pass (run of usable samples), when the dump is complete.
        arrival = (np.arange(n) + 1) * step_s
        if latency_at_pass_end:
            for first, last in _groups(capacity_bits > 0):
                arrival[first : last + 1] = (last + 1) * step_s
        lost_bits = 0.0
        for i in range(n):
            memory += generated_per_step
            if memory_bits and memory > memory_bits:
                # Memory full: the new data are lost (the latency of the rest is then approximate).
                lost_bits += memory - memory_bits
                memory = memory_bits
            take = min(memory, capacity_bits[i])
            if take > 0:
                # First in, first out: the bits downloaded now were generated
                # from (downloaded_total) to (downloaded_total + take) / rate.
                raw.append((take, downloaded_total / generated_bps, (downloaded_total + take) / generated_bps,
                            arrival[i]))
                downloaded_total += take
            memory -= take
            downloaded_bits[i] = take
            memory_end_bits[i] = memory
        result.latency_chunks = _chunks(raw, product_s, product_reference)
        latency = _latency_figures(result.latency_chunks)

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
            **_memory_figures(memory_bits, float(memory_end_bits.max()), lost_bits, generated_bps, n * step_s),
        )
        result.on_board_series = (times, memory_end_bits / 8e6)
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
    antennas: Optional[dict[str, int]] = None,
    setup_s: float = 0.0,
    latency_at_pass_end: bool = False,
    memory_bits: Optional[float] = None,
    product_s: Optional[float] = None,
    product_reference: str = "data",
) -> ConstellationResult:
    """Satellites sharing ground stations whose antennas track one satellite each.

    ``days[k][satellite][station]`` are the samples of day k. Every satellite
    generates ``generated_bps`` into its own memory (empty at the start). At
    each sample each station serves at most one of the satellites for which
    the link closes (margin >= ``min_margin_db``), and each satellite uses at
    most one station: the satellites with more data on board are served
    first, each through its best-margin station still free. A satellite that
    could have downloaded but found its stations busy loses that sample
    (conflict). Latency is first in first out per satellite; the system
    figures put the data of all satellites together. ``antennas`` gives the
    number of antennas of a station (default 1): it serves that many
    satellites at a time. ``setup_s`` is the time an antenna needs to be
    configured for another satellite: after serving one it cannot take a
    different satellite for that long (the same satellite resumes at once).
    ``latency_at_pass_end``: the data downloaded during a pass count as on
    ground at the end of that pass (dump complete), not when each bit comes
    down (the time a product can be processed). ``product_s`` groups each
    satellite's data on ground into products of that many seconds of on-board
    data (see :func:`_chunks`).
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
    n_antennas = np.array([max(1, int((antennas or {}).get(g, 1))) for g in station_names])
    station_usable = usable  # per station, for the contention figures
    # One slot per antenna: the scheduling below gives each slot one satellite.
    slot_station = np.repeat(np.arange(len(station_names)), n_antennas)
    margin, visible, elevation, usable = (a[:, slot_station, :] for a in (margin, visible, elevation, usable))

    n_sat, n_st = len(sat_names), len(station_names)
    step_bits = info_bitrate_bps * step_s
    gen_step = generated_bps * step_s
    memory = np.zeros(n_sat)
    downloaded_total = np.zeros(n_sat)
    served_station = np.full((n_sat, n), -1)
    downloaded = np.zeros((n_sat, n))
    on_board = np.zeros((n_sat, n))
    conflict = np.zeros((n_sat, n), dtype=bool)
    raw: list[list] = [[] for _ in range(n_sat)]  # (bits, generated from, to, on ground), per satellite
    lost = np.zeros(n_sat)
    setup_steps = int(np.ceil(setup_s / step_s)) if setup_s > 0 else 0
    n_slots = len(slot_station)
    last_sat = np.full(n_slots, -1)  # satellite each antenna was last configured for
    free_from = np.zeros(n_slots, dtype=int)  # first sample it may take another satellite
    pending: list[list] = [[] for _ in range(n_sat)]  # chunks of the pass in progress (latency_at_pass_end)

    def _close_pass(s, t_end):
        for take, oldest, newest in pending[s]:
            raw[s].append((take, oldest, newest, t_end))
        pending[s].clear()

    for i in range(n):
        memory += gen_step
        if memory_bits:
            # Memory full: the new data are lost (the latency of the rest is then approximate).
            over = np.maximum(memory - memory_bits, 0.0)
            lost += over
            memory -= over
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
                (
                    g for g in np.flatnonzero(usable[s, :, i])
                    if int(g) not in owner and (last_sat[g] in (-1, s) or i >= free_from[g])
                ),
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
            last_sat[g] = s
            free_from[g] = i + 1 + setup_steps
        for s in range(n_sat):
            if served_station[s, i] < 0:
                if i > 0 and served_station[s, i - 1] >= 0:
                    _close_pass(s, i * step_s)
                continue
            take = min(memory[s], step_bits)
            if take > 0 and generated_bps > 0:
                t_now = (i + 1) * step_s
                oldest = downloaded_total[s] / generated_bps
                newest = (downloaded_total[s] + take) / generated_bps
                if latency_at_pass_end:
                    pending[s].append((take, oldest, newest))
                else:
                    raw[s].append((take, oldest, newest, t_now))
            elif generated_bps <= 0:
                take = step_bits  # capacity only
            memory[s] -= take if generated_bps > 0 else 0.0
            downloaded_total[s] += take
            downloaded[s, i] = take
        on_board[:, i] = memory
    for s in range(n_sat):
        _close_pass(s, n * step_s)
    chunks = [_chunks(raw[s], product_s, product_reference) for s in range(n_sat)]

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
            "lost_mb": float(lost[s]) / 8e6 if with_memory else None,
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
            used = [station_names[slot_station[g]] for g in served_station[s, first : last + 1]]
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

    served_at = np.where(served_station >= 0, slot_station[np.maximum(served_station, 0)], -1)
    for g, name in enumerate(station_names):
        served_here = served_at == g
        result.stations.append({
            "station": name,
            "antennas": int(n_antennas[g]),
            "served_min_per_day": float(np.count_nonzero(served_here)) * step_min / n_days,
            "downloaded_mb_per_day": float(downloaded[served_here].sum()) / 8e6 / n_days,
            # Share of the antenna time in use.
            "busy_pct": float(np.count_nonzero(served_here)) / (n * n_antennas[g]) * 100.0,
            # Minutes with more satellites usable than antennas.
            "contention_min_per_day": float(
                np.count_nonzero(station_usable[:, g, :].sum(axis=0) > n_antennas[g])
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
        result.on_board_series = (times, on_board.sum(axis=0) / 8e6)
        summary.update(
            generated_bps=generated_bps,
            generated_mb_per_day=generated_bps * 86400 / 8e6 * n_sat,
            max_on_board_mb=float(on_board.sum(axis=0).max()) / 8e6,
            final_on_board_mb=float(on_board[:, -1].sum()) / 8e6,
            **_memory_figures(memory_bits, float(on_board.max()), float(lost.sum()), generated_bps, n * step_s),
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
