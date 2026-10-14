# Link Budget Tool

Questo repository contiene un semplice strumento in Python per l'analisi del link budget di collegamenti satellitari. L'applicazione offre un'interfaccia grafica basata su Tkinter che permette di inserire i TLE del satellite, configurare i parametri di trasmissione e visualizzare i risultati.

## Caratteristiche principali
- Predizione delle orbite e delle finestre di contatto tramite Skyfield.
- Calcolo delle perdite di percorso, dell'attenuazione atmosferica ITU e del margine Eb/No.
- Valutazione parallela di downlink e uplink con grafici dedicati per Eb/No e C/No.
- Visualizzazione grafica dell'elevazione e della qualità del link per ogni contatto.
- Funzione di ricalcolo rapido del link budget senza dover ripetere l'analisi delle orbite.
- Finestra principale a schede: **Orbit & Station** (TLE, data, stazione di terra, caricamento dei file), **Downlink Passes** e **Uplink Passes** (parametri e analisi dei passaggi dal TLE), **Static LB · Downlink** e **Static LB · Uplink** (link budget statici).
- Link budget preliminare a elevazione fissa (schede **Static LB**), indipendente dal TLE, nel formato delle tabelle statiche ECSS. Ogni scheda è organizzata in pagine: una barra laterale elenca le pagine degli ingressi e dei risultati (Summary con esito PASS/FAIL e riquadri dei risultati chiave, Link Budget, PFD Check, link budget con tolleranze), i pulsanti Calculate/Load/Save/Export restano sempre visibili in alto, e un riquadro in fondo alla barra mostra l'esito e il margine dell'ultimo calcolo (o che i risultati non sono più aggiornati). In caso di errore di inserimento si apre la pagina del campo sbagliato. Gli ingressi sono divisi nelle sezioni di un link budget statico (stazione e geometria, canale dati e requisiti, trasmettitore, propagazione e limite PFD, ricevitore), con i campi del pannello principale inseriti nella sezione a cui appartengono. Nella sezione del trasmettitore si sceglie esplicitamente la fonte dell'EIRP (valore di EIRP oppure catena Tx: potenza − perdite + guadagno d'antenna) e viene mostrata l'EIRP effettivamente usata; i campi ignorati dal calcolo in quel momento (l'opzione EIRP non scelta, gli axial ratio quando la perdita di polarizzazione è inserita a mano, la larghezza del fascio senza angolo di depuntamento) diventano grigi. Accanto a ogni parametro un'icona ⓘ ne spiega il significato al passaggio del mouse. La perdita ionosferica può essere inserita a mano oppure calcolata dall'indice di scintillazione S4 con il metodo ITU-R P.531 della nota tecnica ESA AWS-TN-ESA-SYS-0052 (S4 riportato alla frequenza del collegamento come f^-1.5, fluttuazione picco-picco 27.5·S4^1.26 dB, fading pari a metà); la finestra mostra la latitudine geomagnetica e la zona ionosferica della stazione (equatoriale, medie latitudini, ovale aurorale, calotta polare) come guida per S4. Per il downlink, la sezione del canale dati può ricavare formattazione e overhead dalla struttura CCSDS TM (frame di 223·I byte con primary header, OCF e FECF facoltativi; ASM a 32 bit per blocco RS(255,223) di 255·I byte; convoluzionale rate 1/2): dato il bit rate codificato del trasmettitore e la profondità di interleaving I, calcola overhead totale, fattore di formattazione (frame / dati utili, togliendo header e trailer del frame e, se si indica la lunghezza dei Source Packet, intestazione, Data Field Header CUC e Packet Error Control di ogni pacchetto) e i bit rate dei transfer frame (riferimento per Eb/N0: tolti convoluzionale, RS e ASM) e informativo (dati utili). Le icone ⓘ accanto ai valori calcolati ne mostrano il dettaglio. La scintillazione troposferica ITU-R P.618 si attiva o esclude nella finestra, e i risultati avvisano quando è fuori dal campo di validità del modello (elevazione sotto 5°, frequenze sotto ~4 GHz) o quando gli Other Attenuations rischiano di contare due volte gli effetti atmosferici e i risultati includono: catena Tx (potenza, perdite di circuito e VSWR, guadagno effettivo, EIRP), perdite di propagazione (spazio libero, ITU-R, ionosfera, polarizzazione calcolabile dagli axial ratio), perdita di puntamento dell'antenna ricevente, S/No, degradazioni, Eb/N0 riferito al bit rate con formattazione, margini statistici (MEAN − 3σ e Worst Case RSS) da tolleranze favorevoli/avverse TRI/UNI/GAU, e il controllo della densità di flusso (PFD) rispetto al limite regolamentare, con la spiegazione del calcolo. Passando il mouse su una riga se ne vede la definizione; i margini sono segnati ✓/✗ rispetto ai requisiti. I risultati si esportano in CSV oppure, con **Export PDF...**, in un report A4 impaginato (sintesi con esito PASS/FAIL dei requisiti, parametri di ingresso e derivati, tabella del link budget, controllo PFD, tolleranze, metodo e ipotesi). I valori condivisi con la scheda dei passaggi della stessa direzione (frequenza, EIRP, G/T, bit rate, roll-off, efficienza spettrale, ecc.) compaiono con l'etichetta in blu, come caselle modificabili collegate: cambiarli in una scheda li cambia anche nell'altra. Gli ingressi si possono salvare/caricare con **Save...** / **Load...** (vedi sotto).

## Requisiti
- Python 3.8 o superiore.
- Dipendenze: `numpy`, `pandas`, `matplotlib`, `skyfield`, `itur`, `astropy`, `scipy`, `reportlab` (per il report PDF), `tkinter` (incluso in Python). È possibile installarle con:
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
3. All'avvio l'applicazione carica automaticamente, se presenti nella cartella del progetto (o accanto all'eseguibile), i file `tle.txt` e `parameters.json`, oltre a `ground_stations.txt` -- proprio come le stazioni di terra. Puoi comunque sovrascriverli in qualsiasi momento con i pulsanti **Load TLE from file**, **Load Ground Stations** e **Import Link Budget Parameters**, oppure modificando i campi a mano.
4. Nella scheda **Downlink Passes** premi **Start/Refresh Analysis** per calcolare le finestre di contatto, quindi seleziona una finestra dalla lista per vedere grafici e tabelle del link budget.

### Formato dei file di ingresso
- **Ground stations (`ground_stations.txt`)**: ogni riga non vuota deve avere 4 campi separati da virgola: `nome, latitudine_deg, longitudine_deg, altitudine_m`. Le righe che iniziano con `#` sono commenti. Esempio:
  ```text
  # name, lat [deg], lon [deg], alt [m]
  Darmstadt, 49.8700, 8.6500, 144
  Lannion, 48.7333, -3.4542, 31
  ```
- **File TLE (`tle.txt`)**: devono contenere almeno due righe consecutive che iniziano con `1 ` e `2 ` (la prima riga del nome è facoltativa). Il pulsante **Load TLE from file** accetta file con tre righe (nome + linee 1-2) o due righe (solo linee 1-2). Se un file `tle.txt` è presente nella cartella del progetto (o accanto all'eseguibile), viene caricato automaticamente all'avvio -- il repository ne include uno di esempio (TLE della ISS) da sostituire con quello del proprio satellite.
- **Parametri di simulazione (`parameters.json`)**: è possibile caricare tutti i campi numerici della GUI da un JSON. Le chiavi supportate sono `eirp_sat_dbw`, `eirp_gs_dbw`, `frequency_ghz`, `uplink_frequency_ghz`, `c_io_dbhz`, `gt_gs_dbk`, `gt_sat_dbk`, `antenna_diameter_m`, `link_availability_pct`, `other_attenuations_db`, `uplink_other_attenuations_db`, `bitrate_mbps`, `rolloff`, `demod_loss_db`, `demod_loss_sat_db`, `overhead`, `spectral_efficiency_bpshz`, `uplink_bitrate_mbps`, `uplink_rolloff`, `uplink_overhead`, `uplink_spectral_efficiency_bpshz`. Il tasso di pioggia per R001 viene stimato automaticamente dalla libreria ITU-R e non deve essere specificato. Se un file `parameters.json` è presente nella cartella del progetto (o accanto all'eseguibile), viene caricato automaticamente all'avvio, come per il file TLE.
  Il file può contenere anche le sezioni facoltative `fixed_elevation_downlink` e `fixed_elevation_uplink` con gli ingressi dei due link budget a elevazione fissa (stazione di terra, elevazione, quota, requisiti, fonte dell'EIRP e catena Tx, perdite, ricevitore, limite PFD e tolleranze). La stazione di queste sezioni è propria delle schede statiche: non cambia quella selezionata nella scheda Orbit & Station, che viene usata solo se la sezione non ne indica una. Non serve scriverle a mano: nelle schede Static LB il pulsante **Save...** le scrive nel file scelto (proposto: il `parameters.json` in uso), insieme ai valori condivisi (etichette blu), mantenendo tutte le altre chiavi; **Load...** le rilegge. Se salvate in `parameters.json`, vengono caricate automaticamente all'avvio insieme agli altri parametri, e anche quando si importa un file di parametri. Il `parameters.json` del repository contiene già la sezione di downlink con i valori della tabella di riferimento.
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
## Creare un eseguibile
Procedura completa su Windows, dalla cartella del progetto aggiornata a `main`. L'ordine conta: la chiave pubblica della licenza deve esistere **prima** di compilare, altrimenti l'eseguibile nasce senza controllo della licenza (vedi [Licenza a scadenza dell'eseguibile](#licenza-a-scadenza-delleseguibile)).

1. **Prepara l'ambiente** (inclusi PyInstaller e `pytest`, richiesto da `astropy.tests.runner` usato durante il packaging):
   ```cmd
   git pull
   pip install -r requirements.txt
   pip install pyinstaller pytest
   ```
2. **Solo la prima volta, se vuoi il controllo della licenza:** crea la coppia di chiavi e metti la chiave pubblica nel repository, così sarà inclusa in tutte le build successive:
   ```cmd
   python license_tool.py keygen --private-key C:\secure\link_budget_license_private_key.pem
   git add assets\license_public_key.pem
   git commit -m "Add licence public key"
   git push
   ```
   Se `assets\license_public_key.pem` esiste già, salta questo passo: rifare le chiavi invalida le licenze già emesse.
3. **Compila** con lo spec già configurato (include `assets\` con logo, icona e chiave pubblica):
   ```cmd
   pyinstaller gui.spec
   ```
   Lo spec crea la cartella `dist\Link-Budget-Tool\`, con `Link-Budget-Tool.exe` e le sue librerie (in `_internal\`). È la modalità "cartella": l'eseguibile parte in pochi secondi perché non deve decomprimere nulla, e appena lanciato mostra una schermata con il logo. Si distribuisce l'intera cartella (ad esempio come zip). In `gui.spec` si può impostare `ONEFILE = True` per avere un unico `.exe`, che però a ogni avvio decomprime alcune centinaia di MB in una cartella temporanea e può impiegare minuti a partire, soprattutto con l'antivirus attivo.
4. **Emetti la licenza** per il cliente (vedi sotto per i dettagli):
   ```cmd
   python license_tool.py issue --private-key C:\secure\link_budget_license_private_key.pem --customer "Nome cliente" --expires 2027-12-31 -o license.key
   ```
5. **Prepara la cartella da consegnare:** in `dist\Link-Budget-Tool\`, accanto a `Link-Budget-Tool.exe`, aggiungi:
   - `license.key`
   - `ground_stations.txt`, `tle.txt`, `parameters.json`

   Poi comprimi l'intera cartella `Link-Budget-Tool` in uno zip.
6. **Prova prima di consegnare:**
   - senza `license.key` nella cartella, `Link-Budget-Tool.exe` deve mostrare "No licence found…" e chiudersi;
   - con `license.key`, in alto a destra deve comparire "Licensed to … · valid until …".

Se la build fallisce o l'eseguibile non parte, conserva il messaggio d'errore completo (log di PyInstaller o traceback) per correggere `gui.spec`.

Il file `__main__.py` è predisposto per essere usato come entry point anche quando il pacchetto è incorporato nell'eseguibile, quindi non servono modifiche aggiuntive.

### File delle ground stations, TLE e parametri
- **Non** vengono inclusi nel bundle PyInstaller: i file `ground_stations.txt`, `tle.txt` e `parameters.json` devono stare nella stessa cartella dell'eseguibile (e dei file Excel che usi insieme al programma).
- All'avvio l'applicazione legge quei file e popola automaticamente il menu delle stazioni di terra, i campi TLE e i parametri di link budget; puoi modificarli o sostituirli senza ricostruire l'eseguibile. Un file mancante o non valido viene semplicemente ignorato (i campi restano vuoti, senza bloccare l'avvio).
- In alternativa puoi specificare un percorso diverso impostando le variabili d'ambiente `GROUND_STATIONS_FILE`, `TLE_FILE` o `PARAMETERS_FILE` prima di avviare il programma.
- Dalla GUI puoi comunque sovrascrivere qualsiasi cosa in qualunque momento con i pulsanti **Load Ground Stations**, **Load TLE from file** e **Import Link Budget Parameters** nella sezione "Observation Settings"; i campi vengono aggiornati immediatamente.

### Licenza a scadenza dell'eseguibile
L'eseguibile può essere limitato con un file di licenza che scade in una data stabilita (`licensing.py`, lo stesso sistema dell'HKTM CCSDS Signal Generator). La licenza è un piccolo file `license.key`, da mettere nella stessa cartella dell'eseguibile, con il nome del cliente e la data di scadenza, firmato con una chiave privata che ha solo chi emette le licenze. L'eseguibile contiene solo la chiave pubblica: può verificare una licenza ma non crearne una, e modificare il file (ad esempio la data) la invalida. Senza una licenza valida il programma si ferma all'avvio con una spiegazione; la licenza valida è indicata nell'intestazione della finestra, in arancione negli ultimi 30 giorni.

1. **Una volta sola:** crea la coppia di chiavi (richiede `cryptography`, già in `requirements.txt`):
   ```cmd
   python license_tool.py keygen --private-key C:\secure\link_budget_license_private_key.pem
   ```
   - La **chiave privata** va dove indichi: tienila al sicuro e fuori dal repository. Chi la possiede può emettere licenze; se va persa servono una nuova coppia di chiavi e un nuovo eseguibile.
   - La **chiave pubblica** va in `assets/license_public_key.pem`. Va inclusa nel repository: ogni eseguibile costruito da quel momento verifica la licenza. Gli eseguibili costruiti senza questo file non hanno il controllo.
2. **Per ogni cliente:** emetti una licenza e inviala insieme all'eseguibile:
   ```cmd
   python license_tool.py issue --private-key C:\secure\link_budget_license_private_key.pem --customer "Nome cliente" --expires 2027-12-31 -o license.key
   ```
   La data di scadenza è l'ultimo giorno valido. Per prolungarla basta emettere un nuovo file, senza ricostruire l'eseguibile.
3. **Verifica facoltativa:** `python license_tool.py verify license.key`

Limiti: è un deterrente, non una protezione dalla copia. L'ultima data vista viene ricordata nel profilo utente, così riportare indietro l'orologio viene rilevato, ma un utente determinato può comunque decompilare un eseguibile Python. Avviando il programma dai sorgenti il controllo non si applica mai.

### Note sugli avvisi di PyInstaller
- Lo spec `gui.spec` include esplicitamente `astropy.tests.runner` e `pytest` per evitare l'errore `ModuleNotFoundError: No module named 'astropy.tests.runner'` durante l'esecuzione del binario.
- Per eliminare l'avviso `WARNING: Hidden import "scipy.special._cdflib" not found!`, lo spec esclude il modulo opzionale `_cdflib` e forza l'inclusione degli altri moduli compilati di `scipy.special`.

## License
This project is released under the [MIT License](LICENSE).
