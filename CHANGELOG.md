# Changelog

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
- Pass analysis moved out of the window into `link_core.py`, used by the window, the statistics and the new command line tool `lbtool.py`.
- Regression tests on a frozen AWS case, validation kit against StatPro (`validation/`), tests on GitHub Actions at every change.

## 1.0.0
- First delivered version: pass analysis, static link budgets (ECSS layout, tolerances, PFD, PDF report), licence check, Windows executable.
