# AWS downlink, Svalbard, 3° elevation: tool vs SatMaster Pro

Static LB · Downlink of `missions/AWS/parameters.json` (600 km, 3°, slant range
2515.5 km, 1.707 GHz, 99.99 %) against a SatMaster Pro NON-GEO link budget
(15 October 2026). SatMaster models a transparent transponder: the uplink was
made ideal (uplink C/N0 101 dB-Hz, 27 dB above the downlink) and the
transponder used at saturation with no back-off.

## SatMaster inputs (downlink)

| Field | Value | Origin in the tool |
|---|---|---|
| Elevation / range | 3° / 2515.5 km | Static LB geometry |
| Antenna gain | +30.7 dBi | G/T 9.9 + 10 log 120 K |
| LNB noise temperature | 94.3 K (1.22 dB) | Tsys 120 K − sky 25.7 K at 5° |
| Antenna noise | 40.4 K | clear sky at 3° |
| Antenna mispoint | 0.01 dB | 0.1° on a 3.97° beam |
| Other path losses | 2.93 dB | ionospheric 2.59 + polarisation 0.09 + radome 0.25 |
| Satellite EIRP | 5.01 dBW | Tx chain 4.51 dBW + pattern 0.5 dBi at 65.9° off nadir |
| Info rate / FEC | 1.556 Mbps / 0.436 | transfer frame rate; 3.57 Mbps coded |
| Modulation | QPSK, 1 + roll-off = 2 | 3.57 MHz occupied |
| Implementation loss | 1.5 dB | demodulator 1.0 + modulation 0.5 |
| C/IM, OBO | 200 dB, 0 dB | single carrier, no intermodulation |

## Results

| | Tool clear | SatMaster clear | Tool rain (p = 0.01 %) | SatMaster rain down |
|---|---|---|---|---|
| Free space loss [dB] | 165.11 | 165.10 | 165.11 | 165.10 |
| Gas [dB] | 0.65 | 0.68 | 0.65 | 0.72 |
| Cloud + rain [dB] | 0.00 | 0.00 | 0.06 | 0.03 |
| Tropospheric scintillation [dB] | — | — | 1.86 (if on) | 2.08 |
| System noise [K] | 134.66 | 137.16 | 137.63 | 141.57 |
| C/N0 [dB-Hz] | **74.31** | **74.21** | **72.35** (scintillation on) | **71.95** |

Clear sky: 0.10 dB, from the system noise (SatMaster adds about 2.7 K) and the
gas absorption (SatMaster: water vapour at 5 % of the time). Rain: 0.40 dB,
from the scintillation (SatMaster at p = 0.0051 %, its share of the overall
availability), the system noise and the gas.

The comparison showed that the tool added sky noise from the tropospheric
scintillation; fixed: only gas, cloud and rain radiate noise (ITU-R P.618 §3).
