# Changelog

## Unreleased
- Statistics: on-board memory per satellite (On-board Memory [GB], lbtool.py --memory-gb): share used, autonomy without contacts, data lost when full; EPS-Sterna: 4 GB (Mission Data partition).

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
