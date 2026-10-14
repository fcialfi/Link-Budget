#!/usr/bin/env python3
"""Link Budget Tool from the command line: passes and statistics without the window.

Uses the missions of satellites.json (parameters.json, antenna pattern,
tle.txt or custom orbit of each mission folder) and the same calculations as
the window (link_core).

Examples:
    python lbtool.py passes --mission AWS --station Svalbard --date 2026-10-14 -o passes.csv
    python lbtool.py statistics --mission AWS --stations Svalbard,Kangerlussuaq \\
        --start 2026-10-14 --days 30 --instruments 100 --hktm 4 -o stats.csv
    python lbtool.py missions
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import date, datetime, timezone

import calculations
import link_core
import missions
import pass_statistics
import tle_source


def _catalogue() -> missions.Catalogue:
    path = calculations.resolve_optional_data_file(missions.CATALOGUE_FILE, env_var="SATELLITES_FILE")
    if not path:
        raise SystemExit("satellites.json not found.")
    return missions.load_catalogue(path)


def _mission(name: str, day: date, plane: str | None, download: bool = False):
    """(satellite, parameters or None, TLE lines) of a mission; loads its antenna pattern."""
    catalogue = _catalogue()
    sat = catalogue.find(name)
    if sat is None:
        names = ", ".join(s.name for s in catalogue.satellites)
        raise SystemExit(f"Unknown mission '{name}'. Missions: {names}")
    params_path = catalogue.parameters_path(sat)
    params = None
    if params_path:
        with open(params_path, encoding="utf-8") as f:
            params = json.load(f)
    pattern = catalogue.pattern_path(sat)
    if pattern:
        calculations.load_antenna_pattern(pattern)
    else:
        calculations.set_isotropic_antenna_pattern()
    if sat.custom_orbit:
        orbit = (params or {}).get("orbit") or {}
        planes = orbit.get("ltan") or ["22:20"]
        planes = planes if isinstance(planes, list) else [planes]
        tle = tle_source.synthetic_tle(
            float(orbit.get("altitude_km", 628.0)),
            tle_source.ltan_hours(plane or planes[0]),
            link_core.utc_midnight(day),
            orbit.get("inclination_deg"),
            name=sat.name.upper(),
            eccentricity=float(orbit.get("eccentricity") or 0.0),
            arg_perigee_deg=float(orbit.get("arg_perigee_deg") or 0.0),
        )
        lines = (tle.line1, tle.line2)
    else:
        path = catalogue.mission_file(sat, missions.TLE_FILE)
        if download and sat.norad_id is not None:
            tle = tle_source.fetch_tle(sat.norad_id)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                f.write(tle.text())
            print(f"TLE downloaded from {tle.source}, saved in {path}", file=sys.stderr)
        try:
            with open(path, encoding="utf-8") as f:
                rows = [line.strip() for line in f if line.strip()]
        except OSError:
            raise SystemExit(f"No TLE for {sat.name}: {path} (use --download, or download it from the window).") from None
        lines = tuple(rows[-2:])
    return sat, params, lines


def _link(params, args):
    if params is None:
        print("No link budget parameters for this mission: contact analysis only.", file=sys.stderr)
        return None
    return link_core.link_inputs_from_parameters(params)


def _stations(names: str, params) -> list[link_core.Station]:
    default_d = float((params or {}).get("antenna_diameter_m") or 3.0)
    result = []
    for name in [n.strip() for n in names.split(",") if n.strip()]:
        if name not in calculations.GROUND_STATIONS:
            raise SystemExit(f"Unknown station '{name}'. Stations: {', '.join(calculations.GROUND_STATIONS)}")
        result.append(link_core.station_from_catalogue(name, default_d))
    return result


def cmd_missions(_args) -> None:
    for sat in _catalogue().satellites:
        print(sat.label)


def cmd_passes(args) -> None:
    day = date.fromisoformat(args.date)
    _sat, params, tle = _mission(args.mission, day, args.plane, args.download)
    link = _link(params, args)
    (station,) = _stations(args.station, params)
    analysis = link_core.PassAnalysis(tle, station, link, args.min_elevation)
    frame = analysis.day(link_core.utc_midnight(day))
    contacts = analysis.contacts(frame)
    if link:
        print(analysis.atmosphere_summary())
    for k, (start, end, direction) in enumerate(contacts, 1):
        in_pass = (frame["Time (UTC)"] >= start) & (frame["Time (UTC)"] <= end)
        text = (f"Contact {k:2d}: {start:%H:%M:%S} - {end:%H:%M:%S} ({direction})  "
                f"max El {frame.loc[in_pass, 'Elevation (°)'].max():5.1f}°")
        if "Eb/No (dB)" in frame:
            text += f"  min Eb/No {frame.loc[in_pass & frame['Visible'], 'Eb/No (dB)'].min():6.2f} dB"
        if "PFD Margin (dB)" in frame:
            text += f"  min PFD margin {frame.loc[in_pass & frame['Visible'], 'PFD Margin (dB)'].min():+5.1f} dB"
        print(text)
    if args.output:
        frame[frame["Elevation (°)"] > 0].to_csv(args.output, index=False)
        print(f"Samples above the horizon written to {args.output}")


def cmd_statistics(args) -> None:
    start = date.fromisoformat(args.start)
    _sat, params, tle = _mission(args.mission, start, args.plane, args.download)
    link = _link(params, args)
    stations = _stations(args.stations, params)
    budget_kwargs = {}
    required_ebno = None
    if link:
        section = params.get(calculations.FIXED_LINK_BUDGET_SECTIONS["downlink"]) or {}
        required_ebno = args.required_ebno if args.required_ebno is not None else section.get("required_ebno_db")
        margin = args.required_margin if args.required_margin is not None else section.get("required_margin_db")
        if required_ebno is None:
            raise SystemExit("Give --required-ebno (no Required Eb/No in the mission parameters).")
        rates = [r for r in (args.instruments, args.hktm, args.other) if r]
        budget_kwargs = dict(
            info_bitrate_bps=link.info_bitrate_bps,
            min_margin_db=float(margin or 0.0),
            generated_bps=sum(rates) * 1e3 if rates else None,
        )
    analyses = [link_core.PassAnalysis(tle, st, link, args.min_elevation) for st in stations]

    def day_samples(day_start):
        per_station = {}
        for analysis in analyses:
            frame = analysis.day(day_start, budget="visible")
            margin_db = frame["Eb/No (dB)"].to_numpy(dtype=float) - float(required_ebno) if link else None
            per_station[analysis.station.name] = pass_statistics.DaySamples(
                day_start, frame["Elevation (°)"].to_numpy(dtype=float), frame["Visible"].to_numpy(dtype=bool),
                margin_db,
            )
        return pass_statistics.combine_stations(per_station)

    def progress(done, total):
        print(f"\r{done}/{total} days", end="", file=sys.stderr, flush=True)

    result = pass_statistics.compute(day_samples, start, args.days, link_core.STEP_S, progress=progress,
                                     **budget_kwargs)
    print(file=sys.stderr)
    for key, value in result.summary.items():
        print(f"{key}: {value:.3f}" if isinstance(value, float) else f"{key}: {value}")
    if args.output:
        with open(args.output, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Date", "Contacts", "Contact [min]", "Usable [min]", "Capacity [MB]",
                             "Downloaded [MB]", "Max Gap [min]", "On Board [MB]"])
            for d in result.days:
                writer.writerow([d.date, d.contacts, d.contact_min, d.usable_min, d.capacity_mb, d.downloaded_mb,
                                 d.max_gap_min, d.on_board_mb])
            writer.writerow([])
            writer.writerow(["Start (UTC)", "End (UTC)", "Stations", "Max El [deg]", "Contact [min]",
                             "Usable [min]", "Capacity [MB]", "Downloaded [MB]", "Min Margin [dB]"])
            for p in result.passes:
                writer.writerow([p.start.isoformat(), p.end.isoformat(), p.stations, p.max_elevation_deg,
                                 p.contact_min, p.usable_min, p.capacity_mb, p.downloaded_mb, p.min_margin_db])
        print(f"Per-day and per-pass tables written to {args.output}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Link Budget Tool without the window.")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("missions", help="list the missions of satellites.json").set_defaults(func=cmd_missions)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--mission", required=True, help="mission name in satellites.json, e.g. AWS")
    common.add_argument("--min-elevation", type=float, default=calculations.MIN_ELEVATION_DEG,
                        help="minimum elevation of a contact [deg] (default 5)")
    common.add_argument("--plane", help="LTAN hh:mm of the orbital plane (custom-orbit missions)")
    common.add_argument("-o", "--output", help="CSV file to write")
    common.add_argument("--download", action="store_true",
                        help="download the latest TLE (SatNOGS / N2YO) into the mission folder first")

    passes = sub.add_parser("passes", parents=[common], help="passes of one day over one station")
    passes.add_argument("--station", required=True)
    passes.add_argument("--date", default=datetime.now(timezone.utc).date().isoformat(), help="YYYY-MM-DD (UTC)")
    passes.set_defaults(func=cmd_passes)

    stats = sub.add_parser("statistics", parents=[common], help="statistics over several days and stations")
    stats.add_argument("--stations", required=True, help="comma-separated, e.g. Svalbard,Kangerlussuaq")
    stats.add_argument("--start", required=True, help="first day, YYYY-MM-DD (UTC)")
    stats.add_argument("--days", type=int, default=7)
    stats.add_argument("--required-ebno", type=float, help="dB (default: mission parameters)")
    stats.add_argument("--required-margin", type=float, help="dB (default: mission parameters)")
    stats.add_argument("--instruments", type=float, help="instrument data rate on the downlink [kbps]")
    stats.add_argument("--hktm", type=float, help="HKTM rate on the downlink [kbps]")
    stats.add_argument("--other", type=float, help="other data rate on the downlink [kbps]")
    stats.set_defaults(func=cmd_statistics)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
