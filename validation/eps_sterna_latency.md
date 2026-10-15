# EPS-Sterna data latency: tool vs reference table

Constellation of `missions/EPS-Sterna/parameters.json` (3 planes × 2 satellites,
phasing 60°, 595.6 km), 64 kbps generated on board per satellite, downlink
1.542 Mbps information rate, minimum elevation 5°, required margin 3 dB, one
antenna per site, 3 days from 2026-10-15. Latency: time from generation on
board to arrival on ground of each bit (first in first out). Reference: table
of stored mission data (SMD) latency provided by the user (SVL3 = Svalbard).

| Network | Reference 50 / 80 / 90 % | Tool 50 / 80 / 90 % [min] |
|---|---|---|
| **Svalbard + Punta Arenas** | **33 / 61 / 77** | **29.7 / 59.6 / 73.5** |
| Svalbard + McMurdo | 24 / 38 / 43 | 18.4 / 33.0 / 38.3 |
| Svalbard + Troll | 27 / 45 / 55 | 21.3 / 39.7 / 53.7 |
| Svalbard + O'Higgins | 29 / 49 / 68 | 23.1 / 43.2 / 61.3 |
| Svalbard only | 47 / 75 / 84 | 52.5 / 85.1 / 139.0 |

McMurdo, Troll and O'Higgins were used for this check only (not in
`ground_stations.txt`): McMurdo −77.8391, 166.6671; Troll −72.0117, 2.5350;
O'Higgins −63.3211, −57.9008.

With two stations the tool is 3-6 min earlier than the reference in every
case, with the same shape: consistent with the reference counting the data as
available at the end of the pass (dump complete), about half a pass later on
average. With Svalbard only the tool is later (90 %: 139 vs 84 min): the
short low-elevation passes usable with a 3 dB margin above 5° cannot empty
one orbit of data and the memory queues for a few orbits; a lower minimum
elevation or no required margin in the reference would explain it.

Open: the reference's definition of arrival (per bit or end of dump), its
minimum elevation, margin and downlink bit rate.

With 104 kbps (an earlier assumption) Svalbard + Punta Arenas gave
32.5 / 63.4 / 75.9 min.

## Minimum pass duration

Contacts shorter than a minimum duration left out (Min Pass Duration), same
cases (50 / 80 / 90 % [min]):

| Network | Reference | 0 min | 2 min | 3 min | 4 min | 5 min |
|---|---|---|---|---|---|---|
| **Svalbard + Punta Arenas** | **33 / 61 / 77** | 29.7 / 59.6 / 73.5 | 30.5 / 61.3 / 74.5 | 31.4 / 62.7 / 75.5 | **32.5 / 64.5 / 77.1** | 34.8 / 68.2 / 80.7 |
| Svalbard + McMurdo | 24 / 38 / 43 | 18.4 / 33.0 / 38.3 | 18.8 / 33.9 / 40.3 | 19.8 / 35.6 / 54.7 | 20.8 / 37.8 / 79.7 | 22.2 / 51.2 / 126.3 |
| Svalbard + Troll | 27 / 45 / 55 | 21.3 / 39.7 / 53.7 | 21.5 / 40.8 / 57.3 | 22.2 / 42.5 / 60.9 | 22.9 / 44.5 / 63.4 | 23.9 / 47.0 / 66.1 |
| Svalbard + O'Higgins | 29 / 49 / 68 | 23.1 / 43.2 / 61.3 | 23.4 / 44.4 / 63.5 | 24.0 / 46.2 / 65.7 | 25.6 / 51.6 / 69.2 | 27.7 / 57.1 / 72.4 |
| Svalbard only | 47 / 75 / 84 | 52.5 / 85.1 / 139.0 | 53.4 / 87.1 / 151.7 | 54.4 / 98.0 / 174.4 | 56.2 / 126.7 / 207.2 | 60.2 / 171.4 / 262.8 |

Leaving out the contacts shorter than 3-4 min brings Svalbard + Punta Arenas
within 1-3 min of the reference on all three percentiles (4 min: 32.5 / 64.5 /
77.1), and Svalbard + O'Higgins and the 80 % of Svalbard + Troll close to it;
the reference probably schedules only passes of a few minutes or more.

With the minimum SDA pass duration of EPS-Sterna, 4.7 min (in the mission
parameters; with the 30 s step a contact needs 10 samples, i.e. it acts as
5 min): Svalbard + Punta Arenas 34.8 / 68.2 / 80.7 min, against 33 / 61 / 77.

## Sensitivity (Svalbard + Punta Arenas, 64 kbps)

| Minimum elevation | Required margin | Min pass 0 min | Min pass 4.7 min |
|---|---|---|---|
| 5° | 3 dB or 0 dB | 29.7 / 59.6 / 73.5 | 34.8 / 68.2 / 80.7 |
| 3° | 3 dB or 0 dB | 25.8 / 53.1 / 69.9 | 28.1 / 58.3 / 72.6 |
| 0° | 3 dB | 23.5 / 49.1 / 66.7 | 24.2 / 51.5 / 68.4 |
| 0° | 0 dB | 22.6 / 47.6 / 65.5 | 23.2 / 49.7 / 67.2 |

- A 10 s step instead of 30 s changes nothing (4.7 min: 34.5 / 67.9 / 80.3).
- The required margin does not matter above 3°: the link closes with more
  than 3 dB in every usable sample.
- The minimum elevation is the dominant assumption: with the 4.7 min minimum
  pass the reference (33 / 61 / 77) lies between 5° and 3°, as with a
  minimum elevation of about 4° or a horizon mask lower than 5° in parts of
  the sky.

Adopted for EPS-Sterna: minimum elevation 3°, minimum pass 4.7 min, 64 kbps:
Svalbard + Punta Arenas 28.1 / 58.3 / 72.6 min (reference 33 / 61 / 77).

## Table 30 (reference) vs tool, all percentiles

Tool with the adopted EPS-Sterna settings: minimum elevation 3°, minimum pass
4.7 min, 64 kbps, 7 days from 2026-10-15 (Inuvik for this check only:
68.3176, −133.5498). Latency [min] at 50 / 60 / 70 / 80 / 90 / 95 / 99 / 100 %.

| Network | Reference | Tool |
|---|---|---|
| SVL | 47 56 65 75 84 89 92 96 | 46 56 66 76 88 126 173 191 |
| SVL+MCM | 24 29 33 38 43 45 48 52 | 16 21 26 31 36 40 77 93 |
| SVL+OHI | 29 35 41 49 68 81 91 95 | 22 28 34 42 62 74 84 90 |
| SVL+TRO | 27 32 38 45 55 74 89 94 | 20 25 31 39 54 70 82 86 |
| **SVL+PUN** | **33 40 49 61 77 85 92 96** | **28 37 46 59 73 80 86 92** |
| SVL+INU | 47 56 65 75 84 89 95 103 | 36 46 55 64 73 78 84 87 |

- Svalbard only: 50-80 % within 1 min, 90 % within 4 min.
- Two stations: the tool is 2-7 min earlier on 50-90 %, the same offset in
  every network, consistent with the reference counting the data at the end
  of the dump (about half a pass later).
- SVL+INU: in the reference Inuvik adds nothing to Svalbard (same figures as
  SVL); in the tool it shortens the latency, since every contact downloads.
  The reference probably downloads once per orbit (one dump, at the first
  station of the orbit), which Inuvik, seeing the same orbits as Svalbard,
  cannot improve.
- Tails (95-100 %): the tool is longer with Svalbard only and SVL+MCM, where
  after an orbit without a contact the next pass cannot download two orbits
  of data at 1.542 Mbps and the memory queues; the reference behaves as if
  every dump emptied the memory.

## Data counted at the end of the dump

Same settings, with the data downloaded during a pass counted as on ground at
the end of that pass (Latency at end of pass), 7 days:

| Network | Reference | Tool, end of pass |
|---|---|---|
| SVL | 47 56 65 75 84 89 92 96 | 54 64 75 85 95 130 178 197 |
| SVL+MCM | 24 29 33 38 43 45 48 52 | 25 30 35 40 45 48 82 99 |
| SVL+OHI | 29 35 41 49 68 81 91 95 | 30 36 43 52 70 83 93 97 |
| SVL+TRO | 27 32 38 45 55 74 89 94 | 28 34 40 48 64 80 92 96 |
| **SVL+PUN** | **33 40 49 61 77 85 92 96** | **36 45 55 67 82 89 95 99** |
| SVL+INU | 47 56 65 75 84 89 95 103 | 48 58 67 77 87 92 96 103 |

Counting the data at the end of the dump reproduces the reference within
1-3 min for SVL+MCM (to 95 %), SVL+OHI (all percentiles), SVL+TRO (to 80 %)
and SVL+INU (all percentiles, 100 % identical), and explains why Inuvik adds
nothing to Svalbard in the reference. The reference therefore counts the
data at the end of the dump, the convention for product timeliness; adopted
for EPS-Sterna (latency_at_pass_end in the mission parameters). SVL+PUN is
3-6 min above the reference and Svalbard only 7-10 min: these two networks
probably differ in another assumption (minimum elevation or horizon mask of
Svalbard and Punta Arenas). A "one dump per pole passage" rule (no new dump
within 30 min of the previous one) changes nothing and was not kept.

Tool with the adopted settings, Svalbard + Punta Arenas: 80 % of the data on
ground within 67 min (85 % within 75 min).
