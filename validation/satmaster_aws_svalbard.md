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

## Results (final run)

| | Tool | SatMaster | Difference |
|---|---|---|---|
| EIRP per carrier [dBW] | 5.01 | 4.95 | 0.06 (residual SatMaster back-off) |
| Free space loss [dB] | 165.11 | 165.10 | 0.01 |
| Atmosphere, clear [dB] | 0.65 | 0.68 | 0.03 |
| System noise, clear [K] | 134.66 | 137.16 | 0.08 dB |
| **C/N0 clear [dB-Hz]** | **74.31** | **74.16** | **0.15** (0.09 at equal EIRP) |
| **Eb/No clear [dB]** | **12.39** | **12.24** | **0.15** |
| Tropospheric scintillation [dB] | 1.86 | 2.13 | 0.27 |
| System noise, rain [K] | 137.63 | 141.85 | 0.13 dB |
| **C/N0 rain down [dB-Hz]** (scintillation on in the tool) | **72.35** | **71.84** | **0.51** (0.45 at equal EIRP) |

Clear sky: 0.09 dB at equal EIRP, from the system noise (SatMaster adds about
2.7 K) and the gas absorption (SatMaster: water vapour at 5 % of the time).
Rain: 0.45 dB at equal EIRP, mostly the tropospheric scintillation, which at 3°
is an extrapolation of ITU-R P.618 (valid from 5°): the two tools differ in the
time percentage and antenna efficiency (SatMaster 40.8 %, tool 50 %); then the
system noise (0.13 dB) and the gas (0.06 dB).

SatMaster settings needed: C/IM 200 dB (AUTO gives 10.65 dB from the transponder
curve at IBO 0), carriers per transponder AUTO (1 adds a 1.8 dB per-carrier
back-off), carrier spacing factor 2 (at least 1 + roll-off), site altitude in km
(0.01), downlink availability 99.99 (AUTO splits the overall availability).

## 10° elevation

Range 1931.2 km, satellite EIRP 5.07 dBW (pattern 0.56 dBi at 64.2° off nadir),
antenna noise 14.5 K (clear sky at 10°). Uplink HPA power fixed at 2.15 W: with
MIN SatMaster lowers the uplink until the Eb/No target is just met, and the
transparent transponder backs off (5.56 dB here); 2.15 W gives IBO 0.

| | Tool | SatMaster | Difference |
|---|---|---|---|
| EIRP per carrier [dBW] | 5.07 | 5.06 | 0.01 |
| Free space loss [dB] | 162.81 | 162.81 | 0.00 |
| Gas, clear [dB] | 0.19 | 0.20 | 0.01 |
| Tropospheric scintillation [dB] | 0.44 | 0.44 | 0.00 |
| System noise, clear / rain [K] | 108.8 / 109.9 | 111.26 / 112.50 | 0.10 dB |
| **C/N0 clear [dB-Hz]** | **78.06** | **77.94** | **0.12** |
| **C/N0 rain down [dB-Hz]** | **77.57** | **77.44** | **0.13** |
| **Eb/No clear [dB]** | **16.14** | **16.02** | **0.12** |

At 10° the scintillation model is valid and the two tools give the same value;
the remaining 0.1 dB, the same at 3° and 10°, is the 2.5-2.7 K SatMaster adds
to the system noise.

## Conclusion

Static downlink budget validated against SatMaster Pro: within 0.1 dB in clear
sky (3° and 10°), 0.13 dB in rain at 10°, 0.45 dB in rain at 3° where the
tropospheric scintillation is extrapolated below 5°.

The comparison showed that the tool added sky noise from the tropospheric
scintillation; fixed: only gas, cloud and rain radiate noise (ITU-R P.618 §3).
