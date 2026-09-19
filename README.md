# simulazione_completa.py

Script Python per l'esecuzione di simulazioni di dinamica molecolare (MM/MD) 
a partire da un file di configurazione `parameters.toml`.

## Requisiti

È necessario eseguire il programma in un ambiente conda con i seguenti pacchetti installati:
- OpenMM
- OpenFF Toolkit
- OpenFF Interchange
- PACKMOL
- RDKit
- pandas

## Utilizzo

### Calcolo ab initio

Per costruire e avviare una simulazione da zero:

```bash
python /path/to/simulazione_completa.py -f parameters.toml --from_toml True
```

### Ripresa di una simulazione

Per riprendere una simulazione da uno stato salvato:

```bash
python /path/to/simulazione_completa.py -f parameters.toml --from_state True
```

> **Nota:** per riprendere correttamente la simulazione è necessario che i file 
> `final_state.xml` e `interchange.json` siano presenti nella cartella di esecuzione, 
> oppure che il loro percorso sia esplicitato nel file `parameters.toml`.

---

## Output

Al termine della simulazione viene creata una cartella `run_YYYYMMDD_hhmmss` 
nella directory del file `parameters.toml`. Una copia degli stessi file viene 
salvata nella cartella `Last_Run`, che viene sovrascritta ad ogni ripresa della simulazione.

| File | Descrizione |
|---|---|
| `box.csv` | Step e modulo dei vettori di cella (Lx, Ly, Lz in nm) |
| `data.csv` | Step, energia potenziale, cinetica e totale (kJ/mol), temperatura (K), volume (nm³) e densità (g/mL) |
| `data_full.csv` | Combinazione di `box.csv` e `data.csv` |
| `final_state.xml` | Stato finale della simulazione, usato per riprendere il calcolo |
| `interchange.json` | Oggetto `openMM.Interchange` corrispondente allo stato finale |
| `parameters.toml` | Copia del file di configurazione usato |
| `run.log` | Log completo dell'esecuzione |
| `topology.pdb` | Topologia iniziale costruita da PACKMOL |
| `trajectory.pdb` | Traiettoria della simulazione |

---

## File di configurazione `parameters.toml`

Il file è diviso in cinque sezioni.

### `[molecules_paths]`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `HBD_path` | str | Percorso del file `HBD.sdf` |
| `HBA_path` | str | Percorso del file `HBA.sdf` |
| `counterion_path` | str | Percorso del file `CI.sdf` |
| `cellobiose_topology` | str | Percorso del file PDB della topologia del cellobiosio |

### `[molecules_copies]`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `HBD_copies` | int | Numero di molecole HBD |
| `HBA_copies` | int | Numero di molecole HBA |
| `counterion_copies` | int | Numero di molecole di controione |
| `water_copies` | int | Numero di molecole di H₂O |

### `[openff_methods]`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `partial_charges_method` | str | Metodo per le cariche parziali: `am1bcc`, `am1bccelf10`, `am1-mulliken`, `mmff94`, `gasteiger` |
| `box_padding` | float | Padding attorno alla cella di cellobiosio (nm) |
| `simulation_steps` | int | Numero di passi di simulazione (1 step = 1 fs) |
| `packmol_random_seed` | int | Seed per PACKMOL |

### `[simulation_parameters]`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `simulation_temperature` | float | Temperatura del termostato (K) |

### `[resume_simulation]`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `interchange_path` | str | Percorso del file `interchange.json` |
| `state_path` | str | Percorso del file `final_state.xml` |
| `simulation_steps` | int | Numero di passi da eseguire alla ripresa (1 step = 1 fs) |

---

## Esempio di `parameters.toml`

```toml
[molecules_paths]
HBD_path = "/path/to/HBD.sdf"
HBA_path = "/path/to/HBA.sdf"
counterion_path = "/path/to/CI.sdf"
cellobiose_topology = "/path/to/cellobiose.pdb"

[molecules_copies]
HBD_copies = 136
HBA_copies = 68
counterion_copies = 68
water_copies = 228

[openff_methods]
partial_charges_method = "am1bcc"
box_padding = 20.0
simulation_steps = 500000
packmol_random_seed = 11472

[simulation_parameters]
simulation_temperature = 350.0

[resume_simulation]
interchange_path = "./Last_Run/interchange.json"
state_path = "./Last_Run/final_state.xml"
simulation_steps = 5000
```
