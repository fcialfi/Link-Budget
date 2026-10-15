# Changelog

## Unreleased
- Statistics of a whole constellation (custom orbit: planes × satellites per plane, phasing): shared stations tracking one satellite at a time, conflicts, per-satellite and per-station figures; also in lbtool.py (--constellation).
- Statistics: use of each station of the network.
- Statistics: latency requirement (share of the data on ground within a time).
- Orbit & Passes: Constellation row with editable planes (LTAN list), satellites per plane and phasing; choice of the satellite analysed (plane, sat #); Show Constellation (table, ground tracks, satellites in their planes). lbtool.py: --sat.
- Statistics plots redrawn: contact time, data per day, gaps (or time per satellite) and data latency with the requirement.
- Sun outage: Sun noise in the ground antenna beam along the passes (Sun position, dish pattern, solar flux), lowering G/T and Eb/No; Sun Separation and Sun Noise columns, passes marked in the contact list and in lbtool.py, included in the statistics.
- Fix: tropospheric scintillation no longer adds sky noise (only gas, cloud and rain absorption radiate noise, ITU-R P.618 §3).
- Validation: AWS downlink at Svalbard against SatMaster Pro (validation/satmaster_aws_svalbard.md): C/N0 within 0.1 dB in clear sky.
- Fix: crash of the statistics on Windows (garbage collection in the worker thread).

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
