# Link Budget Tool

Questo repository contiene un semplice strumento in Python per l'analisi del link budget di collegamenti satellitari. L'applicazione offre un'interfaccia grafica basata su Tkinter che permette di inserire i TLE del satellite, configurare i parametri di trasmissione e visualizzare i risultati.

## Caratteristiche principali
- Predizione delle orbite e delle finestre di contatto tramite Skyfield.
- Calcolo delle perdite di percorso, dell'attenuazione atmosferica ITU e del margine Eb/No.
- Valutazione parallela di downlink e uplink con grafici dedicati per Eb/No e C/No.
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
3. Inserisci manualmente (o carica dai pulsanti **Load TLE from file**, **Load Ground Stations** e **Load Parameters**) le righe TLE, la stazione di terra e gli altri parametri, quindi premi **Start** per calcolare le finestre di contatto.
4. Seleziona una finestra dalla lista per vedere grafici e tabelle del link budget.

### Formato dei file di ingresso
- **Ground stations (`ground_stations.txt`)**: ogni riga non vuota deve avere 4 campi separati da virgola: `nome, latitudine_deg, longitudine_deg, altitudine_m`. Le righe che iniziano con `#` sono commenti. Esempio:
  ```text
  # name, lat [deg], lon [deg], alt [m]
  Darmstadt, 49.8700, 8.6500, 144
  Lannion, 48.7333, -3.4542, 31
  ```
- **File TLE**: devono contenere almeno due righe consecutive che iniziano con `1 ` e `2 ` (la prima riga del nome è facoltativa). Il pulsante **Load TLE from file** accetta file con tre righe (nome + linee 1-2) o due righe (solo linee 1-2).
- **Parametri di simulazione (`parameters.json`)**: è possibile caricare tutti i campi numerici della GUI da un JSON. Le chiavi supportate sono `eirp_sat_dbw`, `eirp_gs_dbw`, `frequency_ghz`, `uplink_frequency_ghz`, `c_io_dbhz`, `gt_gs_dbk`, `gt_sat_dbk`, `antenna_diameter_m`, `link_availability_pct`, `other_attenuations_db`, `uplink_other_attenuations_db`, `bitrate_mbps`, `rolloff`, `demod_loss_db`, `demod_loss_sat_db`, `overhead`, `spectral_efficiency_bpshz`, `uplink_bitrate_mbps`, `uplink_rolloff`, `uplink_overhead`, `uplink_spectral_efficiency_bpshz`. Il tasso di pioggia per R001 viene stimato automaticamente dalla libreria ITU-R e non deve essere specificato.
  Esempio:
  ```json
  {
    "frequency_ghz": 1.707,
    "bitrate_mbps": 3.57,
    "rolloff": 0.45,
    "spectral_efficiency_bpshz": 2.0,
    "uplink_bitrate_mbps": 1.0,
    "uplink_rolloff": 0.35,
    "uplink_overhead": 1.1,
    "uplink_spectral_efficiency_bpshz": 1.5
  }
  ```

## Spunti per personalizzazioni
- Modifica la variabile `GROUND_STATIONS` nel file Python per aggiungere nuove stazioni di terra.
- Adatta la tabella dell'antenna (`ANTENNA_PATTERN_ANGLES` e `ANTENNA_PATTERN_GAINS`) al tuo specifico diagramma.
- Per cambiare la durata dell'analisi o l'intervallo di campionamento, modifica i parametri all'inizio della funzione `run_analysis`.
- Imposta il campo "Spectral Efficiency" nel pannello Baseband per ottenere la larghezza di banda del canale desiderata.
- Usa il pulsante **Export Cesium View** per generare un bundle HTML/CZML e visualizzare la traiettoria e la stazione di terra in 3D su CesiumJS.
- Se hai installato `cesiumpy`, puoi usare **Preview Cesium (CesiumPy)** per aprire subito una preview senza scegliere un percorso di export.

### Export Cesium (CZML + HTML)
È disponibile un export in formato CZML con un viewer HTML pronto all'uso. Il pulsante **Export Cesium View**:

1. Usa i dati di analisi già calcolati per il contatto selezionato.
2. Genera un file `.czml` e un file `.html` nella posizione scelta.
3. Apre automaticamente l'HTML nel browser.

L'export funziona senza dipendenze extra, ma se hai già installato `cesiumpy` verrà utilizzato in modo opzionale. Il viewer HTML usa CesiumJS da CDN, quindi serve una connessione internet quando apri l'HTML.

### Preview Cesium (CesiumPy)
Se `cesiumpy` è disponibile, il pulsante **Preview Cesium (CesiumPy)** crea un bundle temporaneo e apre il viewer automaticamente nel browser. È utile per un'anteprima rapida, mentre l'export serve per condividere i file HTML/CZML o conservarli in una cartella specifica.

> Nota: alcune versioni di CesiumPy richiedono una Google Maps API key. In tal caso, imposta la variabile d'ambiente `CESIUMPY_GOOGLE_API_KEY` prima di avviare il programma (verrà propagata anche alle variabili `GOOGLE_API_KEY` e `GOOGLE_MAPS_API_KEY`).
Se non l'hai impostata, il pulsante **Preview Cesium (CesiumPy)** ti chiederà la key direttamente.

## Creare un eseguibile
Per generare un eseguibile standalone con PyInstaller (ad esempio su Windows):

1. Installa i pacchetti necessari (inclusi PyInstaller e `pytest`, richiesto da `astropy.tests.runner` usato durante il packaging):
   ```bash
   pip install -r requirements.txt
   pip install pyinstaller pytest
   ```
2. Dalla cartella del progetto esegui PyInstaller usando lo spec già configurato:
   ```bash
   pyinstaller gui.spec
   ```
3. Al termine troverai l'eseguibile in `dist/gui/` (su Windows, `dist/gui/gui.exe`).

Il file `__main__.py` è predisposto per essere usato come entry point anche quando il pacchetto è incorporato nell'eseguibile, quindi non servono modifiche aggiuntive.

### File delle ground stations
- **Non** viene incluso nel bundle PyInstaller: il file `ground_stations.txt` deve stare nella stessa cartella dell'eseguibile (e dei file Excel che usi insieme al programma).
- All'avvio l'applicazione legge quel file e popola il menu a tendina delle stazioni di terra; puoi modificarlo o sostituirlo senza ricostruire l'eseguibile.
- In alternativa puoi specificare un percorso diverso impostando la variabile d'ambiente `GROUND_STATIONS_FILE` prima di avviare il programma.
- Dalla GUI puoi selezionare rapidamente un catalogo personalizzato con il pulsante **Load Ground Stations** nella sezione "Observation Settings"; l'elenco viene aggiornato immediatamente.

### Note sugli avvisi di PyInstaller
- Lo spec `gui.spec` include esplicitamente `astropy.tests.runner` e `pytest` per evitare l'errore `ModuleNotFoundError: No module named 'astropy.tests.runner'` durante l'esecuzione del binario.
- Per eliminare l'avviso `WARNING: Hidden import "scipy.special._cdflib" not found!`, lo spec esclude il modulo opzionale `_cdflib` e forza l'inclusione degli altri moduli compilati di `scipy.special`.

## License
This project is released under the [MIT License](LICENSE).
