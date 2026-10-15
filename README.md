# Link Budget Tool

Questo repository contiene un semplice strumento in Python per l'analisi del link budget di collegamenti satellitari. L'applicazione offre un'interfaccia grafica basata su Tkinter che permette di inserire i TLE del satellite, configurare i parametri di trasmissione e visualizzare i risultati.

**Manuale utente:** [docs/Link_Budget_Tool_User_Manual.pdf](docs/Link_Budget_Tool_User_Manual.pdf), in inglese come l'interfaccia. Spiega schede, ingressi, risultati, metodo di calcolo, file di dati e licenza. Per rigenerarlo dopo una modifica: `python docs/manual/take_screenshots.py`, che rifà gli screenshot e richiede uno schermo (su Linux senza schermo va lanciato con `xvfb-run`), poi `python docs/manual/build_manual.py`.

## Caratteristiche principali
- Predizione delle orbite e delle finestre di contatto tramite Skyfield.
- Calcolo delle perdite di percorso, dell'attenuazione atmosferica ITU e del margine Eb/No.
- Valutazione parallela di downlink e uplink con grafici dedicati per Eb/No e C/No.
- Visualizzazione grafica dell'elevazione e della qualità del link per ogni contatto.
- Funzione di ricalcolo rapido del link budget senza dover ripetere l'analisi delle orbite.
- Nome del progetto (campo **Project Name** nella scheda Orbit & Passes, chiave `project_name` di `parameters.json`): compare nel titolo della finestra, nel report PDF (fascia in alto e sopra il titolo), nel nome proposto per il PDF, negli export CSV (in quello dei passaggi insieme al passaggio esportato), e viene salvato con **Save...** delle schede statiche.
- Finestra principale a schede: **Orbit & Passes** (TLE, data, stazione di terra, caricamento dei file, **Start/Refresh Analysis** e lista dei passaggi, comuni a downlink e uplink), **Downlink** e **Uplink** (parametri e risultati per il passaggio selezionato), **Static LB · Downlink** e **Static LB · Uplink** (link budget statici).
- Link budget preliminare a elevazione fissa (schede **Static LB**), indipendente dal TLE, nel formato delle tabelle statiche ECSS. Ogni scheda è organizzata in pagine: una barra laterale elenca le pagine degli ingressi e dei risultati (Summary con esito PASS/FAIL e riquadri dei risultati chiave, Link Budget, PFD Check, link budget con tolleranze), i pulsanti Calculate/Load/Save/Export restano sempre visibili in alto, e un riquadro in fondo alla barra mostra l'esito e il margine dell'ultimo calcolo (o che i risultati non sono più aggiornati). In caso di errore di inserimento si apre la pagina del campo sbagliato. Gli ingressi sono divisi nelle sezioni di un link budget statico (stazione e geometria, canale dati e requisiti, trasmettitore, propagazione e limite PFD, ricevitore), con i campi del pannello principale inseriti nella sezione a cui appartengono. Nella sezione del trasmettitore si sceglie esplicitamente la fonte dell'EIRP (valore di EIRP oppure catena Tx: potenza − perdite + guadagno d'antenna) e viene mostrata l'EIRP effettivamente usata; i campi ignorati dal calcolo in quel momento (l'opzione EIRP non scelta, gli axial ratio quando la perdita di polarizzazione è inserita a mano, la larghezza del fascio senza angolo di depuntamento) diventano grigi. Accanto a ogni parametro un'icona ⓘ ne spiega il significato al passaggio del mouse. La perdita ionosferica può essere inserita a mano oppure calcolata dall'indice di scintillazione S4 con il metodo ITU-R P.531 della nota tecnica ESA AWS-TN-ESA-SYS-0052 (S4 riportato alla frequenza del collegamento come f^-1.5, fluttuazione picco-picco 27.5·S4^1.26 dB, fading pari a metà); la finestra mostra la latitudine geomagnetica e la zona ionosferica della stazione (equatoriale, medie latitudini, ovale aurorale, calotta polare) come guida per S4. Per il downlink, la sezione del canale dati può ricavare formattazione e overhead dalla struttura CCSDS TM (frame di 223·I byte con primary header, OCF e FECF facoltativi; ASM a 32 bit per blocco RS(255,223) di 255·I byte; convoluzionale rate 1/2): dato il bit rate codificato del trasmettitore e la profondità di interleaving I, calcola overhead totale, fattore di formattazione (frame / dati utili, togliendo header e trailer del frame e, se si indica la lunghezza dei Source Packet, intestazione, Data Field Header CUC e Packet Error Control di ogni pacchetto) e i bit rate dei transfer frame (riferimento per Eb/N0: tolti convoluzionale, RS e ASM) e informativo (dati utili). Le icone ⓘ accanto ai valori calcolati ne mostrano il dettaglio. La scintillazione troposferica ITU-R P.618 si attiva o esclude nella finestra, e i risultati avvisano quando è fuori dal campo di validità del modello (elevazione sotto 5°, frequenze sotto ~4 GHz) o quando gli Other Attenuations rischiano di contare due volte gli effetti atmosferici e i risultati includono: catena Tx (potenza, perdite di circuito e VSWR, guadagno effettivo, EIRP), perdite di propagazione (spazio libero, ITU-R, ionosfera, polarizzazione calcolabile dagli axial ratio), perdita di puntamento dell'antenna ricevente, S/No, degradazioni, Eb/N0 riferito al bit rate con formattazione, margini statistici (MEAN − 3σ e Worst Case RSS) da tolleranze favorevoli/avverse TRI/UNI/GAU, e il controllo della densità di flusso (PFD) rispetto al limite regolamentare, con la spiegazione del calcolo. Passando il mouse su una riga se ne vede la definizione; i margini sono segnati ✓/✗ rispetto ai requisiti. I risultati si esportano in CSV oppure, con **Export PDF...**, in un report A4 impaginato (sintesi con esito PASS/FAIL dei requisiti, parametri di ingresso e derivati, tabella del link budget, controllo PFD, tolleranze, metodo e ipotesi). I valori condivisi con la scheda Downlink o Uplink (frequenza, EIRP, G/T, bit rate, roll-off, efficienza spettrale, ecc.) compaiono con l'etichetta in blu, come caselle modificabili collegate: cambiarli in una scheda li cambia anche nell'altra. Gli ingressi si possono salvare/caricare con **Save...** / **Load...** (vedi sotto).

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
3. All'avvio l'applicazione carica automaticamente, se presenti nella cartella del progetto (o accanto all'eseguibile), i file `tle.txt` e `parameters.json`, oltre a `ground_stations.csv` -- proprio come le stazioni di terra. Puoi comunque sovrascriverli in qualsiasi momento con i pulsanti **Load TLE from file**, **Load Ground Stations** e **Import Link Budget Parameters**, oppure modificando i campi a mano.
4. Nella scheda **Orbit & Passes** premi **Start/Refresh Analysis** per calcolare le finestre di contatto, quindi seleziona una finestra dalla lista: le schede **Downlink** e **Uplink** mostrano grafici e tabelle del link budget di quel passaggio.
   - **Perdite comuni con il link budget statico:** le schede Downlink e Uplink, con la casella *Include the losses entered in the Static LB tab* (attiva di default), sottraggono le stesse perdite indipendenti dall'elevazione del link budget statico della stessa direzione: ionosferica, depolarizzazione e multipath (il puntamento dell'antenna del satellite è dato dal diagramma, quello della stazione è nullo). Aggiungono anche la degradazione di modulazione alla perdita del demodulatore e riferiscono Eb/N0 al bit rate dei transfer frame (fattore di formattazione). Il riquadro mostra i valori applicati. Attenuazione atmosferica, EIRP, G/T e le altre perdite restano quelle della scheda dei passaggi.
   - **Diagramma d'antenna del satellite:** è il guadagno assoluto (dBi) dell'antenna del satellite in funzione dell'angolo dal nadir (boresight). Ogni missione ha il suo, nella propria cartella (`missions/<missione>/antenna_pattern.txt`, vedi sotto): non c'è un pattern comune. Il file di AWS (ed EPS-Sterna ne ha una copia) contiene un esempio generico (5,6 dBi al nadir, 0,5 dBi da ±65° a ±66°; oltre la tabella il guadagno è 0 dBi): va sostituito con il diagramma reale dell'antenna. Una missione senza pattern usa un'antenna isotropa (0 dBi), segnalata in arancione. Si può anche caricare un altro file con **Load Antenna Pattern...** nella scheda Orbit & Passes. Si applica sempre, istante per istante:
     - **downlink:** EIRP = catena Tx della scheda Static LB · Downlink (potenza − perdite di circuito − perdita VSWR) + guadagno del diagramma verso la stazione. Se quella scheda usa un valore di EIRP, l'EIRP SAT si intende al picco del diagramma e se ne sottrae il calo verso la stazione;
     - **uplink:** il G/T SAT si intende al picco del diagramma e se ne sottrae il calo;
     - **stazione di terra:** l'antenna insegue il satellite, quindi il suo puntamento non introduce perdite.
     Nella scheda Static LB · Downlink la casella *Antenna gain from the satellite antenna pattern at this elevation* legge il guadagno Tx dallo stesso diagramma, all'angolo dal nadir che corrisponde a elevazione e quota inserite (ad esempio 65,4° a 3° e 628 km). Così statico e passaggi usano lo stesso guadagno, e la Tx Antenna Pointing Loss non serve più.
     Nella scheda Static LB · Uplink la casella *G/T SAT at the antenna boresight: subtract the pattern roll-off at this elevation* fa lo stesso per la ricezione. Il G/T SAT si specifica di solito al boresight (nadir), dove il guadagno è massimo. Verso una stazione a bassa elevazione l'antenna la vede lontano dal nadir, con un guadagno più basso del calo del diagramma (picco − guadagno a quell'angolo): l'opzione sottrae quel calo dal G/T SAT, esattamente come fa l'analisi dei passaggi istante per istante. In quel caso il depuntamento dell'antenna Rx non si usa, perché è già compreso nel diagramma.
   - **Scintillazione troposferica (ITU-R P.618):** è disattivata di default sia nei passaggi (casella *Tropospheric scintillation* nella scheda Downlink) sia nelle schede statiche. Il modello vale da 5° di elevazione e tra ~4 e 20 GHz; in banda L/S a bassa elevazione conta la scintillazione ionosferica, presa dall'indice S4 delle schede statiche.

### Formato dei file di ingresso
- **Ground stations (`ground_stations.csv`)**: CSV con intestazione, si apre e si salva in Excel (anche con il punto e virgola). Colonne: `name, latitude_deg, longitude_deg, altitude_m` (obbligatorie) e, facoltative (celle vuote), `antenna_diameter_m, gt_dbk, gt_elevation_deg, tsys_k, antennas`. Si modifica anche dalla finestra con **Edit Stations...** (tab Orbit & Passes). Il vecchio `ground_stations.csv` senza intestazione è ancora letto. Esempio:
  ```text
  name,latitude_deg,longitude_deg,altitude_m,antenna_diameter_m,gt_dbk,gt_elevation_deg,tsys_k,antennas
  Darmstadt,49.87,8.65,144,,,,,
  Svalbard,78.2232,15.6267,10,3,,,,3
  ```
- **Parametri di missione (`missions/<nome>/parameters.json`)**: si vedono e si modificano come tabella, con unità e descrizioni, con **Mission Parameters...** accanto al menu Satellite.
- **File TLE (`tle.txt`)**: devono contenere almeno due righe consecutive che iniziano con `1 ` e `2 ` (la prima riga del nome è facoltativa). Il pulsante **Load TLE from file** accetta file con tre righe (nome + linee 1-2) o due righe (solo linee 1-2). Se un file `tle.txt` è presente nella cartella del progetto (o accanto all'eseguibile), viene caricato automaticamente all'avvio -- il repository include il TLE di AWS (Arctic Weather Satellite, NORAD 60543). Sotto i campi TLE compaiono il file, l'epoca e l'età del TLE (in arancione se ha più di 7 giorni).
- **Satelliti e missioni**: in cima alla scheda Orbit & Passes il menu **Satellite** elenca i satelliti di `satellites.json` (missioni EUMETSAT, file modificabile: `name`, `norad_id`, `description`, oppure `"custom_orbit": true` per una missione non ancora lanciata). Ogni satellite ha una cartella `missions/<nome>/` con `parameters.json` (parametri di link budget), il pattern d'antenna (facoltativo: `antenna_pattern.txt` oppure un qualsiasi file .txt/.csv con "pattern" nel nome, ad esempio `AWS Antenna Pattern.csv`. Due colonne, angolo dal nadir in gradi e guadagno in dBi, separate da virgola, tab, spazi o punto e virgola, e in questo caso anche con la virgola decimale come scrive Excel in italiano; righe di intestazione e commenti `#` ignorati. Se il file manca o non si legge, il programma non prende il pattern di un'altra missione: usa un'antenna isotropa 0 dBi e lo segnala in arancione sotto i pulsanti del pattern) e `tle.txt` (ultimo TLE scaricato, riscaricato alla selezione se ha più di un giorno; per una custom orbit, l'ultimo TLE sintetico generato). Se `parameters.json` c'è, i parametri vengono caricati e il link budget è attivo; se manca, il programma fa solo l'**analisi dei contatti** (orari, durata ed elevazione massima dei passaggi) e disattiva le schede Downlink, Uplink e Static LB. Per una missione con **custom orbit** (ad esempio EPS-Sterna, di cui AWS è il protoflight) si inseriscono gli elementi medi dell'orbita: altitudine (semiasse maggiore medio meno raggio equatoriale), eccentricità, argomento del perigeo ω, inclinazione (vuota = eliosincrona, calcolata) e LTAN (per una costellazione il menu elenca i piani orbitali). Il programma genera un TLE sintetico con epoca alla data dell'analisi e moto medio scelto in modo che SGP4 abbia il periodo nodale della teoria J2, così gli orari dei passaggi non derivano nel corso della giornata; poi calcola passaggi e link budget come per un satellite lanciato. EPS-Sterna usa la sua orbita di riferimento: 595,6 km, e = 0,001214, ω = 90°, i = 97,794°, piani a 15:30, 19:30 e 23:30 (periodo nodale 96,716 min). L'orbita si salva nella sezione `orbit` del `parameters.json` della missione (`altitude_km`, `eccentricity`, `arg_perigee_deg`, `inclination_deg` o `null`, `ltan` come lista di "hh:mm", il primo è quello usato). Senza `satellites.json` il programma usa, come prima, `tle.txt` e `parameters.json` nella sua cartella.
- **G/T della stazione**: il G/T GS è quello dell'antenna di riferimento (Svalbard, 3 m: 9,9 dB/K), all'elevazione di riferimento e a cielo sereno. Con **G/T Reference Diameter** (3 m) una stazione con antenna diversa (5ª colonna di `ground_stations.csv`, oppure il campo Antenna Diameter GS) ha G/T + 20·log10(D/Dref): stesso feed e ricevitore, 2,4 m = −1,94 dB. Con **System Noise Temperature** (120 K) il G/T segue il rumore di cielo: T_sky = 275 K·(1 − 10^(−A/10)) + 2,7 K·10^(−A/10) dall'attenuazione ITU-R A all'elevazione e con il meteo effettivi (ITU-R P.618), e G/T = G/T_ref − 10·log10((Tsys − T_sky,ref + T_sky)/Tsys); più basso a bassa elevazione e con pioggia. Lo spill-over dal terreno non è modellato. **Radome Loss** (0,25 dB) si sottrae al G/T. Nell'analisi dei passaggi attenuazione atmosferica e G/T sono calcolati a ogni istante, all'elevazione effettiva; nel link budget statico all'elevazione impostata, con le righe di rumore di cielo e di sistema.
- **Elevazione minima e maschera d'orizzonte**: nella scheda Orbit & Passes **Min Elevation** (5° di default) è l'elevazione minima di un contatto. Ogni stazione può avere una maschera d'orizzonte in `horizon_masks/<stazione>.txt` (azimut dal Nord e elevazione minima, lineare tra i punti; vedi `horizon_masks/README.txt`): un campione è in contatto solo sopra entrambe. **Horizon Mask** la mostra in un grafico.
- **PFD in funzione dell'elevazione**: nella scheda Static LB · Downlink **PFD Limit 0–5° El** e, facoltativo, **PFD Limit 25–90° El**, con rampa lineare tra 5° e 25° come nel Regolamento Radio ITU art. 21 (senza il secondo il limite è costante). Il link budget statico usa il limite alla sua elevazione; l'analisi dei passaggi calcola la PFD in 4 kHz a ogni istante (EIRP dal pattern verso la stazione, banda occupata da bit rate, roll-off ed efficienza spettrale) e il margine sul limite: colonne PFD nella tabella e margine minimo di ogni passaggio nella lista dei contatti, con l'avviso PFD EXCEEDED.
- **Statistiche su più giorni e più stazioni** (scheda **Statistics**): per il satellite, l'elevazione minima e le maschere della scheda Orbit & Passes e il link budget della scheda Downlink, somma l'analisi dei passaggi su N giorni e su una rete di stazioni (in contatto se almeno una lo vede; a ogni istante la stazione con il margine migliore; ognuna con il suo diametro e la sua maschera): contatti e minuti di contatto al giorno, minuti utili (margine = Eb/No − Required Eb/No ≥ Required Margin, dalla scheda Static LB · Downlink se i campi sono vuoti), buco massimo e medio tra i contatti, capacità di downlink (bit rate informativo × minuti utili). Con i ratei in kbps di quello che passa sulla banda di downlink (**Instrument Data**, **HKTM**, **Other Data**; l'uplink dei telecomandi non conta) simula la memoria di bordo: volume davvero scaricato, dati a bordo a fine giornata e massimi, latenza dei dati (media e massima, FIFO), e se il downlink tiene il passo. I ratei si leggono e si salvano nel `parameters.json` della missione (`downlink_data_kbps`). Con **Latency Requirement** (min) e **Within** (% dei dati), per esempio 60 min per il 90%, verifica il requisito di latenza: quota dei dati a terra entro il tempo richiesto e tempo entro cui arriva la quota richiesta (chiave `latency_requirement`). Tabelle per giorno e per passaggio (con le stazioni) e uso di ogni stazione, grafico ed export CSV. **Costellazione** (missioni con custom orbit, per esempio EPS-Sterna): con *Whole constellation* simula tutti i satelliti dei piani (lista LTAN) × *Sats/plane*, equispaziati nel piano e sfasati di *Phasing* da un piano al successivo; ogni satellite ha la sua memoria, ogni stazione segue un satellite alla volta e lo tiene finché il link chiude, le stazioni libere vanno ai satelliti con più dati a bordo. Dà conflitti (satellite pronto ma stazioni occupate), volume e latenza di sistema e per satellite, uso delle stazioni. Per un satellite senza parametri solo contatti e buchi.
- **Download TLE**: inserisci il **NORAD ID** e premi **Download TLE** (scheda Orbit & Passes): il programma scarica l'ultimo TLE provando le fonti in ordine finché una risponde: SatNOGS (`db.satnogs.org`, dati Space-Track; conosce solo i satelliti con un trasmettitore radio nel suo database) e N2YO (il TLE della pagina del satellite). Il TLE trovato viene salvato nel `tle.txt` della missione; se nessuna fonte risponde, la riga di stato riporta il motivo per ciascuna. Senza connessione resta il TLE già caricato. Il NORAD ID si salva in `parameters.json` (chiave `norad_id`). Il download usa le impostazioni proxy di sistema.
- **Parametri di simulazione (`parameters.json`)**: è possibile caricare tutti i campi numerici della GUI da un JSON. Le chiavi supportate sono `project_name` (testo), `norad_id` (numero NORAD del satellite, per **Download TLE**), `eirp_sat_dbw`, `eirp_gs_dbw`, `frequency_ghz`, `uplink_frequency_ghz`, `c_io_dbhz`, `min_elevation_deg`, `gt_gs_dbk`, `gt_ref_diameter_m`, `gt_ref_elevation_deg`, `tsys_ref_k`, `radome_loss_db`, `gt_sat_dbk`, `antenna_diameter_m`, `link_availability_pct`, `other_attenuations_db`, `uplink_other_attenuations_db`, `bitrate_mbps`, `rolloff`, `demod_loss_db`, `demod_loss_sat_db`, `overhead`, `spectral_efficiency_bpshz`, `uplink_bitrate_mbps`, `uplink_rolloff`, `uplink_overhead`, `uplink_spectral_efficiency_bpshz`. Un valore `null` lascia il campo vuoto: per `c_io_dbhz` (C/Io) significa nessun interferente, e il link budget usa solo C/No. Il tasso di pioggia per R001 viene stimato automaticamente dalla libreria ITU-R e non deve essere specificato. Se un file `parameters.json` è presente nella cartella del progetto (o accanto all'eseguibile), viene caricato automaticamente all'avvio, come per il file TLE.
  Il file può contenere anche le sezioni facoltative `fixed_elevation_downlink` e `fixed_elevation_uplink` con gli ingressi dei due link budget a elevazione fissa (stazione di terra, elevazione, quota, requisiti, fonte dell'EIRP e catena Tx, perdite, ricevitore, limite PFD e tolleranze). La stazione di queste sezioni è propria delle schede statiche: non cambia quella selezionata nella scheda Orbit & Passes, che viene usata solo se la sezione non ne indica una. Non serve scriverle a mano: nelle schede Static LB il pulsante **Save...** le scrive nel file scelto (proposto: il `parameters.json` in uso), insieme ai valori condivisi (etichette blu), mantenendo tutte le altre chiavi; **Load...** le rilegge. Se salvate in `parameters.json`, vengono caricate automaticamente all'avvio insieme agli altri parametri, e anche quando si importa un file di parametri. Il `parameters.json` del repository contiene già la sezione di downlink con i valori della tabella di riferimento.
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
- Adatta il diagramma d'antenna del satellite predefinito (`ANTENNA_PATTERN_ANGLES` in gradi dal nadir e `ANTENNA_PATTERN_GAINS` in dBi) al tuo, oppure caricalo da file con **Load Antenna Pattern...**.
- Per cambiare la durata dell'analisi o l'intervallo di campionamento, modifica i parametri all'inizio della funzione `run_analysis`.
- Imposta il campo "Spectral Efficiency" nel pannello Baseband per ottenere la larghezza di banda del canale desiderata.
## Riga di comando (`lbtool.py`)
Gli stessi calcoli della finestra, senza finestra, per le missioni di `satellites.json` (dai sorgenti, con il `.venv` attivo):
```cmd
python lbtool.py missions
python lbtool.py passes --mission AWS --station Svalbard --date 2026-10-14 -o passaggi.csv
python lbtool.py statistics --mission AWS --stations Svalbard,Kangerlussuaq --start 2026-10-14 --days 30 --instruments 100 --hktm 4 -o statistiche.csv
python lbtool.py passes --mission EPS-Sterna --plane 19:30 --station Svalbard
```
`--download` scarica prima il TLE aggiornato, `--min-elevation` cambia l'elevazione minima, `--required-ebno` / `--required-margin` sostituiscono i valori della missione, `--latency` / `--latency-pct` danno il requisito di latenza, `--constellation` simula tutta la costellazione. `python lbtool.py statistics -h` elenca tutte le opzioni.

## Struttura del codice e test
- `calculations.py`: modelli (geometria, ITU-R, link budget, G/T, PFD, pattern d'antenna).
- `link_core.py`: analisi dei passaggi di una giornata (visibilità con maschera, link budget, PFD, contatti) e lettura dei parametri di una missione; usata da finestra, statistiche e riga di comando.
- `pass_statistics.py`: statistiche su più giorni e stazioni, memoria di bordo e latenza.
- `tle_source.py` (download e TLE sintetici), `missions.py` (catalogo e cartelle), `static_report.py` (righe e PDF del link budget statico), `gui.py` (finestra).
- Test: `pytest` (unitari, più `tests/test_link_core.py` che confronta un caso AWS congelato in `tests/data/aws_reference`). Girano automaticamente su GitHub a ogni modifica (`.github/workflows/tests.yml`). La validazione con StatPro è in `validation/`.
- Le modifiche di ogni versione sono in `CHANGELOG.md`.

## Creare un eseguibile
**Modo rapido:** fai doppio clic su `build.bat` nella cartella del progetto.
- **Prima volta:** crea il `.venv` e ci installa i requirements e PyInstaller; le volte successive aggiunge solo quello che manca.
- **Build:** cancella `build` e `dist`, compila e copia accanto all'eseguibile `ground_stations.csv`, `satellites.json` e la cartella `missions`.
- **Zip (facoltativo):** di default non viene creato. Per consegnare il programma lancia `build.bat zip` dal prompt dei comandi: crea anche `dist\Link-Budget-Tool_<versione>.zip`, da consegnare dopo averci aggiunto il `license.key` del cliente.

I passi sotto descrivono lo stesso procedimento a mano, compresa la creazione delle chiavi di licenza (una volta sola, prima della prima build con licenza).

**Versione:** il numero di versione sta solo in `version.py` (`APP_VERSION`). Compare nel titolo e nell'intestazione della finestra, nella schermata di avvio, nel piè di pagina del report PDF e nel nome dello zip (`build.bat zip`). Per una nuova consegna:
1. aggiorna `APP_VERSION`;
2. fai commit;
3. crea il tag con lo stesso numero (`git tag v1.0.0`, poi `git push --tags`).

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
   - `ground_stations.csv`, `satellites.json` e la cartella `missions` (con le sottocartelle)

   Poi comprimi l'intera cartella `Link-Budget-Tool` in uno zip.
6. **Prova prima di consegnare:**
   - senza `license.key` nella cartella, `Link-Budget-Tool.exe` deve mostrare "No licence found…" e chiudersi;
   - con `license.key`, in alto a destra deve comparire "Licensed to … · valid until …".

Se la build fallisce o l'eseguibile non parte, conserva il messaggio d'errore completo (log di PyInstaller o traceback) per correggere `gui.spec`.

Il file `__main__.py` è predisposto per essere usato come entry point anche quando il pacchetto è incorporato nell'eseguibile, quindi non servono modifiche aggiuntive.

### File delle ground stations, TLE e parametri
- **Non** vengono inclusi nel bundle PyInstaller: i file `ground_stations.csv`, `satellites.json` e la cartella `missions` devono stare nella stessa cartella dell'eseguibile (e dei file Excel che usi insieme al programma).
- All'avvio l'applicazione legge quei file e popola automaticamente il menu delle stazioni di terra, i campi TLE e i parametri di link budget; puoi modificarli o sostituirli senza ricostruire l'eseguibile. Un file mancante o non valido viene semplicemente ignorato (i campi restano vuoti, senza bloccare l'avvio).
- In alternativa puoi specificare un percorso diverso impostando le variabili d'ambiente `GROUND_STATIONS_FILE`, `TLE_FILE`, `PARAMETERS_FILE` o `ANTENNA_PATTERN_FILE` prima di avviare il programma.
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

## Licenza
Software proprietario, tutti i diritti riservati: vedi [LICENSE](LICENSE). L'uso dell'eseguibile richiede una licenza cliente valida (`license.key`).
