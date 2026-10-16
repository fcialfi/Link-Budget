# Changelog

## Unreleased
- Statistics: ground products. Product Length [min] (product_length_min, lbtool.py --product-length): the ground processes the data in products of that many minutes of on-board data (fixed blocks from 00:00 UTC), ready when the last of their data is on ground; only complete products count. Product Latency From (product_latency_reference, --product-reference): each datum's sensing, or the start or end of the product's sensing. EPS-Sterna: 30 min.
- User manual Issue 1.5: new chapters 10 (Parameter and calculation reference) and 11 (Statistics reference): every input with its parameters.json key, unit, value used when empty and the formula through which it reaches the results; every result of the pass analysis, the static link budget and the statistics with its formula and the function that computes it; worked examples on an invented mission and station (G/T versus elevation in clear sky and rain, with the effect of the reference elevation; static S-band downlink; memory and latency of the statistics) computed by the program when the manual is built. The manual text no longer names AWS. Formulas typeset as mathematics (TeX notation drawn as vector outlines with the STIX fonts), in displayed equation blocks with notes and in the formula columns of the tables; chapter 9 rewritten with displayed equations. List of symbols (section 10.2): every symbol of the formulas with its meaning, unit and the section where it is defined. The parameters.json table lists every key with name, unit and meaning; chapter references updated.
- Constellation: Extra sats (orbit.extra_satellites, "LTAN@deg", e.g. 19:30@120) adds single satellites to the regular planes, in the statistics, Show Constellation, 3D View, the satellite selector and lbtool.py --plane/--sat.
- parameters.json in a standard form: every known parameter, in the same order for every mission, null when unused (e.g. AWS shows where to enter the instrument data rate); every save writes this form. A null sun_noise means on.
- Ground stations in ground_stations.csv, a CSV with a header that opens in Excel (semicolon separator accepted; the old ground_stations.txt is still read); Edit Stations... edits the catalogue as a table.
- Mission Parameters...: the parameters.json of the selected mission as a table grouped by subject, with units and descriptions, edited and saved from the window (mission_schema.py).
- Statistics: "Contacts shown" also for a single satellite (one day in detail), with a note saying why it is greyed out (no statistics yet, or one day only).
- Statistics, whole constellation: the contact timeline covers the whole period (one row per satellite, the days on the time axis) instead of the first day only; "Contacts shown" picks one day in detail.
- Fix: selecting another mission clears every parameter its parameters.json does not give (link, statistics, station network, static link budgets back to their defaults), and all of them for a mission without parameters.json, instead of keeping the values of the mission before.
- Fix: the statistics results are cleared when another satellite (mission, or plane / sat # for a single-satellite run) is selected, instead of staying on screen.
- missions/EPS-Sterna/README.md: data sheet of the EPS-Sterna values (orbit, link, on-board data, memory, ground segment, requirement) with their sources and the latest results.
- Statistics: Save to Mission and Reload buttons: the statistics inputs are saved in and reloaded from the parameters.json of the selected mission.

## 1.4.0

### Statistics
- On-board memory per satellite (On-board Memory [GB], lbtool.py --memory-gb): largest share used, autonomy without contacts, data lost when full.
- Real-time telemetry during the passes (Real-time TM [kbps], lbtool.py --realtime), taken from the downlink bit rate available for the stored data.

### EPS-Sterna
- L-band downlink data: payload 55.076 kbps (instrument data rate budget), stored HKTM 3.105 kbps (388.1 B/s), events and system logs 1.216 kbps (882 kB per orbit), real-time HKTM 4.616 kbps (577 B/s) during the passes.
- On-board memory 4 GB (Mission Data partition).

## 1.3.0

### Statistics
- Minimum pass duration (Min Pass Duration [min], lbtool.py --min-pass): contacts shorter than it are not used, per satellite and station, also across midnight.
- Latency at end of pass (lbtool.py --latency-at-pass-end): the data of a pass count as on ground when the dump is complete, the convention for product timeliness.
- Plots without on-board data rates: revisit gaps of each satellite and contact time by station.
- Fix: the on-board data rates and the latency requirement of a mission are cleared when another mission without them is selected.

### EPS-Sterna
- Mission parameters: 64 kbps generated on board, minimum SDA pass 4.7 min, minimum elevation 3°, latency at end of pass. AWS: minimum elevation 5°.
- lbtool.py takes the minimum elevation from the mission parameters unless --min-elevation is given.

### Validation
- EPS-Sterna data latency against the reference SMD latency table (validation/eps_sterna_latency.md): with the data counted at the end of the dump, within 1-3 min for SVL+MCM, SVL+OHI, SVL+TRO and SVL+INU; sensitivity to time step, minimum elevation, margin and minimum pass.

## 1.2.0

### Constellation
- Statistics of a whole constellation (custom orbit: planes × satellites per plane, phasing): stations tracking one satellite per antenna, conflicts, per-satellite and per-station figures; also in lbtool.py (--constellation).
- Orbit & Passes: Constellation row with editable planes (LTAN list), satellites per plane and phasing; choice of the satellite analysed (plane, sat #); Show Constellation (table, ground tracks, satellites in their planes) and 3D View with a time slider. lbtool.py: --sat.
- EPS-Sterna: phasing 60° between the planes (the 6 satellites 60° apart).
- Antenna setup time between satellites (Antenna Setup [min], lbtool.py --setup; 0 by default, 7 min is the requirement's worst case).

### Ground stations
- Own measured G/T (with its elevation and Tsys) and number of antennas per station in ground_stations.txt; one antenna per site is used, the others are backups, unless Use backup antennas (--backup-antennas) is checked.
- Svalbard with 3 antennas; Punta Arenas (Ruta 9 Norte km 28, 2 antennas) added.
- Statistics network per mission (statistics_stations): Svalbard + Punta Arenas for the global data of AWS and EPS-Sterna.

### Link budget
- Sun outage: Sun noise in the ground antenna beam along the passes (Sun position, dish pattern, solar flux), lowering G/T and Eb/No; Sun Separation and Sun Noise columns, passes marked in the contact list and in lbtool.py, included in the statistics.
- Doppler rate (Hz/s) along the passes, from the range acceleration: table column, Doppler plot, largest shift and rate per pass in the contact list and in lbtool.py.
- Fix: tropospheric scintillation no longer adds sky noise (only gas, cloud and rain absorption radiate noise, ITU-R P.618 §3).

### Statistics
- Latency requirement (share of the data on ground within a time) and use of each station of the network.
- New plots: contact timeline (per day, or per satellite for a constellation) coloured by station, data on board over time, data latency with the requirement.
- PDF report of the statistics (Export PDF…): verdict, inputs, summary, plots, tables per station, satellite, day and pass.
- Fix: crash of the statistics on Windows (garbage collection in the worker thread).

### Validation
- AWS downlink at Svalbard against SatMaster Pro (validation/satmaster_aws_svalbard.md): C/N0 within 0.12 dB in clear sky at 3° and 10°, 0.13 dB in rain at 10°.

## 1.1.0

### Missions and orbits
- Satellite menu (`satellites.json`, EUMETSAT missions) with one folder per mission in `missions/`: parameters, antenna pattern, last TLE. A satellite without parameters gets the contact analysis only.
- TLE download by NORAD ID (SatNOGS, then N2YO), with epoch and age of the TLE.
- Custom orbit for missions not launched yet, from mean elements (altitude, eccentricity, argument of perigee, inclination or sun-synchronous, LTAN of each plane); the synthetic TLE keeps the J2 nodal period. EPS-Sterna with its reference orbit.
- Antenna pattern of each mission only from its folder (any `*pattern*.txt/csv`, Excel CSV accepted); isotropic 0 dBi when missing.

### Link budget
- Atmospheric attenuation at the elevation of every sample of the passes.
- Ground station G/T versus elevation and weather (sky noise, ITU-R P.618), radome loss, scaling with the antenna diameter (5th column of `ground_stations.txt`).
- PFD limit versus elevation (ITU RR Art. 21 shape), checked along the passes; 1700–1710 MHz limits of Table 21-4 for AWS and EPS-Sterna.
- Horizon mask per station (`horizon_masks/`) and editable minimum elevation.

### Statistics
- Statistics tab: contacts, usable minutes, capacity and gaps over several days and a network of ground stations; on-board memory model (instrument, HKTM and other data rates) with downloaded volume, memory occupation and data latency.

### Usability
- Pass table with short headings, fitted columns and column groups.
- Author and About box (F1); project name only in the window title.
- `build.bat`: zip only with `build.bat zip`, clear message when `dist` is in use.
- User manual (PDF) with screenshots.

### Quality
- Proprietary licence (LICENSE), replacing the MIT licence of the repository.
- Pass analysis moved out of the window into `link_core.py`, used by the window, the statistics and the new command line tool `lbtool.py`.
- Regression tests on a frozen AWS case, validation kit against StatPro (`validation/`), tests on GitHub Actions at every change.

## 1.0.0
- First delivered version: pass analysis, static link budgets (ECSS layout, tolerances, PFD, PDF report), licence check, Windows executable.
