# EPS-Sterna – mission data sheet

Values used by the tool for EPS-Sterna (`parameters.json` in this folder), with
their source. They are loaded when EPS-Sterna is selected; change them in the
window and save (Static LB Save, or Save to Mission in the Statistics tab).

## Orbit and constellation

| Value | Source |
|---|---|
| 595.6 km mean altitude, e 0.001214, ω 90°, i 97.794° | EPS-Sterna reference orbit table |
| Planes at LTAN 15:30, 19:30, 23:30 | reference orbit table |
| 2 satellites per plane, phasing 60° (the six satellites 60° apart) | user; with 0° the planes reach the pole together and conflict |

## Link (L-band downlink, 1.707 GHz)

As AWS, the protoflight (same link parameters, see `missions/AWS`): Tx chain
3.3 W, circuit loss 0.5 dB, VSWR 1.5, antenna pattern of this folder, coded
bit rate 3.57 Mbps (information 1.542 Mbps), QPSK, roll-off 1.0, required
Eb/No 2.6 dB, required margin 3 dB, PFD limits of RR Table 21-4.
Ground station: G/T 9.9 dB/K at 5° clear sky for the 3 m Svalbard antenna,
Tsys 120 K, radome 0.25 dB.

## On-board data on the L-band downlink

| Flow | Rate | Source |
|---|---|---|
| Payload (instrument) | 55.076 kbps | instrument data rate budget: 62006 bits per scan, one scan every 1.126 s (8.5 km at 7.55 km/s); includes the instrument telemetry |
| Stored HKTM | 3.105 kbps | 388.1 B/s |
| Events and system logs | 1.216 kbps | 882 kB per orbit (96.7 min) |
| Real-time HKTM | 4.616 kbps, during the passes | 577 B/s; also on the L-band downlink |

Payload and HKTM go on different virtual channels of the same downlink; the
tool treats the stored data as one first-in first-out queue.

## On-board memory

64 GB per band, partitions of 4 GB (1024 blocks of 4 MB; a block is 128 pages
of 32 kB). The tool uses 4 GB, the Mission Data partition.

## Ground segment and operations

| Value | Source |
|---|---|
| Global data network: Svalbard + Punta Arenas | user |
| One antenna per site used (Svalbard has 3, Punta Arenas 2; the others are backups) | user |
| Punta Arenas site: Ruta 9 Norte km 28 (−52.93, −70.85) | site address |
| Minimum elevation 3° | chosen against the reference SMD latency table |
| Minimum SDA pass duration 4.7 min | operations |
| Antenna setup between satellites: 0 (requirement: station configured for a pass in less than 7 min, from scratch) | requirement; 7 min is the worst case |
| Latency counted at the end of the pass (dump complete) | matches the reference SMD latency table |

## Requirements

- Timeliness: 80 % of the L1 products disseminated within 75 min of sensing
  (end-to-end: on-board storage, downlink, ground transfer, L1 processing,
  dissemination).

## Results with these values (tool 1.4.0, 2026-10-15, 3 to 7 days)

| Figure | Value |
|---|---|
| Data latency on ground, Svalbard + Punta Arenas, 50 / 80 / 90 % | 35.8 / 66.9 / 81.6 min |
| Reference SMD latency table (Table 30), SVL+PUN | 33 / 61 / 77 min |
| Time left for the ground chain at 80 % (75 − 67) | about 8 min (reference: about 14) |
| Conflicts between satellites | 0 min/day |
| On-board memory used | about 1 % of 4 GB; full after about 6 days without contacts |

Validation: `validation/eps_sterna_latency.md` (latency) and
`validation/satmaster_aws_svalbard.md` (static link budget).
