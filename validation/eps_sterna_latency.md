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
