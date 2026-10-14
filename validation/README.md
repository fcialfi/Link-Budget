# Validazione con StatPro

Confronto del link budget statico di downlink di AWS (Svalbard) con StatPro.

## File

- `AWS_statpro_cases.csv`: i casi, già calcolati dal nostro programma (colonne *Ours*), con le colonne *StatPro* vuote da riempire.
- `make_statpro_cases.py`: rigenera il file con i parametri attuali di `missions/AWS/parameters.json` (serve uno schermo; su Linux `xvfb-run`). Va rilanciato se i parametri cambiano.
- `compare_statpro.py`: confronta le due colonne e elenca le differenze oltre 0,1 dB (`--tolerance` per cambiare la soglia).

## Casi

| Caso | Elevazione | Note |
|---|---|---|
| 1 | 3° | caso di progetto, con la tolleranza sul guadagno Tx (±1,5 dB, TRI) |
| 2 | 10° | |
| 3 | 90° | zenit |

Altitudine 600 km. Colonne *Clear Sky* (solo gas) e *Rain Faded* (pioggia, nubi alla disponibilità del 99,99%).

## Ingressi da mettere in StatPro

| Parametro | Valore |
|---|---|
| Frequenza | 1707 MHz |
| Stazione | Svalbard (78,2232° N, 15,6267° E, 10 m), antenna 3 m |
| Altitudine del satellite | 600 km |
| EIRP | quella della riga *EIRP [dBW]* di ogni caso (Tx 3,3 W − 0,5 dB circuito − perdita VSWR 1,5 + guadagno dal pattern a quell'elevazione) |
| G/T stazione | 9,9 dB/K costante (nel confronto il G/T in funzione dell'elevazione è spento) |
| Disponibilità | 99,99% |
| Perdita ionosferica | riga *Ionospheric Loss* (da S4 = 0,3 a 1575,42 MHz) |
| Perdita di polarizzazione | 0,09 dB |
| Puntamento Rx | depuntamento 0,1°, fascio 3,97° |
| Perdita del demodulatore | 1,0 dB; degrado di modulazione 0,5 dB |
| Bit rate per Eb/No | riga *Transfer Frame Bit Rate [dBHz]* (61,92 dBHz: CCSDS, interleaving 5) |
| Eb/No richiesto | 2,6 dB; margine richiesto 3 dB |
| PFD | banda occupata 3570 kHz (3,57 Mbps, roll-off 1, 2 b/s/Hz), limite −154 dB(W/m²) in 4 kHz |

## Procedura

1. Apri `AWS_statpro_cases.csv` in Excel e scrivi i valori di StatPro nelle colonne *StatPro - Clear Sky* e *StatPro - Rain Faded*, riga per riga (lascia vuote le righe che StatPro non ha).
2. Salva come CSV (va bene anche con il punto e virgola di Excel in italiano).
3. `python validation/compare_statpro.py`
4. Le differenze oltre 0,1 dB vanno spiegate (modelli diversi di atmosfera, ionosfera, polarizzazione…) o corrette.
