# Link Budget Tool

Questo repository contiene un semplice strumento in Python per l'analisi del link budget di collegamenti satellitari. L'applicazione offre un'interfaccia grafica basata su Tkinter che permette di inserire i TLE del satellite, configurare i parametri di trasmissione e visualizzare i risultati.

## Caratteristiche principali
- Predizione delle orbite e delle finestre di contatto tramite Skyfield.
- Calcolo delle perdite di percorso, dell'attenuazione atmosferica ITU e del margine Eb/No.
- Visualizzazione grafica dell'elevazione e della qualità del link per ogni contatto.
- Funzione di ricalcolo rapido del link budget senza dover ripetere l'analisi delle orbite.

## Requisiti
- Python 3.8 o superiore.
- Dipendenze: `numpy`, `pandas`, `matplotlib`, `skyfield`, `itur`, `astropy`, `scipy`, `tkinter` (incluso in Python). È possibile installarle con:
  ```bash
  pip install -r requirements.txt
  ```

## Utilizzo
1. Clona il repository e spostati nella cartella:
   ```bash
   git clone <url_del_repository>
   cd link_budget
   ```
2. Avvia l'applicazione con:
   ```bash
   python -m link_budget
   ```
3. Inserisci le due righe TLE, scegli la stazione di terra e gli altri parametri, quindi premi **Start** per calcolare le finestre di contatto.
4. Seleziona una finestra dalla lista per vedere grafici e tabelle del link budget.

## Spunti per personalizzazioni
- Modifica la variabile `GROUND_STATIONS` nel file Python per aggiungere nuove stazioni di terra.
- Adatta la tabella dell'antenna (`ANTENNA_PATTERN_ANGLES` e `ANTENNA_PATTERN_GAINS`) al tuo specifico diagramma.
- Per cambiare la durata dell'analisi o l'intervallo di campionamento, modifica i parametri all'inizio della funzione `run_analysis`.
- Imposta il campo "Spectral Efficiency" nel pannello Baseband per ottenere la larghezza di banda del canale desiderata.

## Creare un eseguibile
Per generare un eseguibile standalone con PyInstaller (ad esempio su Windows):

1. Installa i pacchetti necessari (incluso PyInstaller):
   ```bash
   pip install -r requirements.txt
   pip install pyinstaller
   ```
2. Dalla cartella del progetto esegui PyInstaller usando lo spec già configurato:
   ```bash
   pyinstaller gui.spec
   ```
3. Al termine troverai l'eseguibile in `dist/gui/` (su Windows, `dist/gui/gui.exe`).

Il file `__main__.py` è predisposto per essere usato come entry point anche quando il pacchetto è incorporato nell'eseguibile, quindi non servono modifiche aggiuntive.

### Note sugli avvisi di PyInstaller
- Lo spec `gui.spec` esclude i moduli di test di Astropy e il pacchetto `pytest` per evitare l'avviso `ModuleNotFoundError: No module named 'pytest'`.
- Per eliminare l'avviso `WARNING: Hidden import "scipy.special._cdflib" not found!`, lo spec esclude anche il modulo opzionale `_cdflib` e forza l'inclusione degli altri moduli compilati di `scipy.special`.

## License
This project is released under the [MIT License](LICENSE).

