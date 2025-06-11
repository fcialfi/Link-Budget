# Link Budget Tool

Questo repository contiene un semplice strumento in Python per l'analisi del link budget di collegamenti satellitari. L'applicazione offre un'interfaccia grafica basata su Tkinter che permette di inserire i TLE del satellite, configurare i parametri di trasmissione e visualizzare i risultati.

## Caratteristiche principali
- Predizione delle orbite e delle finestre di contatto tramite Skyfield.
- Calcolo delle perdite di percorso, dell'attenuazione atmosferica ITU e del margine Eb/No.
- Visualizzazione grafica dell'elevazione e della qualità del link per ogni contatto.
- Funzione di ricalcolo rapido del link budget senza dover ripetere l'analisi delle orbite.

## Requisiti
- Python 3.8 o superiore.
- Dipendenze: `numpy`, `pandas`, `matplotlib`, `skyfield`, `itur`, `astropy`, `scipy`, `tkinter` (incluso in Python). È possibile installare i pacchetti mancanti con:
  ```bash
  pip install numpy pandas matplotlib skyfield itur astropy scipy
  ```

## Utilizzo
1. Clona il repository e spostati nella cartella:
   ```bash
   git clone <url_del_repository>
   cd Link-Budget
   ```
2. Avvia l'applicazione con:
   ```bash
   python 'Link Budget Tool'
   ```
3. Inserisci le due righe TLE, scegli la stazione di terra e gli altri parametri, quindi premi **Start** per calcolare le finestre di contatto.
4. Seleziona una finestra dalla lista per vedere grafici e tabelle del link budget.

## Spunti per personalizzazioni
- Modifica la variabile `GROUND_STATIONS` nel file Python per aggiungere nuove stazioni di terra.
- Adatta la tabella dell'antenna (`ANTENNA_PATTERN_ANGLES` e `ANTENNA_PATTERN_GAINS`) al tuo specifico diagramma.
- Per cambiare la durata dell'analisi o l'intervallo di campionamento, modifica i parametri all'inizio della funzione `run_analysis`.
