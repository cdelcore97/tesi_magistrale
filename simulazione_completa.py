'''
Versione del 10/3/2026
Cristian Del Core
Novità:
- Cambio barostato in MonteCarlo ANISOTROPIC. Permette di bloccare una o più dimensioni della cella unitaria (bloccato asse z)
- Sistema con stesse dimensioni fisiche, ma con più molecole di solvente all'interno per evitare implosioni
- Misura delle dimensioni della box e output in linea (non nel .csv)
- Stampa finale di un paio di grafici con pandas + matplotlib.pyplot
- Stampa in un unico file del csv con dimensioni della simulazione
- Semplificazione input/output. Input in un file parameters.toml

Colli di bottiglia noti:
- Creazione delle cariche sulla cellulosa (~28 minuti)
- Creazione dell'interchange **molto variabile** (20 minuti o 20 secondi. Indaga)

Da fare:
- Output nel file .csv di pressione e dimensioni della simulazione
- Parametrizzazione del sistema (più semplice da impostare)
'''

import time
import datetime
import logging
from pathlib import Path
import pandas as pd
import sys
import os
import tomllib      #Importa parametri da un file parameters.toml. Contiene path di molecole e altri parametri
import argparse
import shutil
import mdtraj
import nglview
import numpy as np
import openmm
import openmm.app
import openmm.unit
from openff.toolkit import ForceField, Molecule, unit, Topology
from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields
from rich.pretty import pprint
from openff.toolkit.utils import get_data_file_path

from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box

#per calcolare energie
from openff.interchange.drivers import get_openmm_energies
from openff.interchange.drivers.all import get_summary_data
from rdkit import Chem

# WORK IN PROGRESS
## Classe reporter per calcolare la pressione interna del sistema mediante teorema del viriale

""" class ReporterEsteso(openmm.app.StateDataReporter):
    #aggiunge al csv standard anche Lx, Ly, Lz (dimensioni del box in nm) e la pressione interna (bar)

    def __init__ (self, file, reportInterval, **kwargs):
        kwargs['volume'] = True
        super().__init__(file, reportInterval, **kwargs)
        self._needsPositions = False
        self._needsVelocities = False
        self._needEnergy = True 
        self._needsForces = True #Per la pressione

    def _constructHeaders(self):
        headers = super()._constructHeaders()
        headers += ['Lx (nm)', 'Ly (nm)', 'Lz (nm)', 'Pressione (bar)']
        return headers
    
    def _constructReportValues(self, simulation, state):
        values = super()._constructReportValues(simulation, state)

        #dimensioni box
        box = state.getPeriodicBoxVectors(asNumpy = True).value_in_unit(openmm.unit.nanometer)
        Lx = box[0,0]
        Ly = box[1,1]
        Lz = box[2,2]

        # Pressione interna tramite il tensore del viriale
        # P = (2*KE + W) / (3*V)  dove W è il viriale = sum(r_i · f_i)
        state_full = simulation.context.getState(
            getPositions=True,
            getForces=True,
            getEnergy=True,
            enforcePeriodicBox=True,
        )
        
        #pressione interna
        positions = state_full.getPositions(asNumpy=True).value_in_unit(openmm.unit.nanometer)
        forces = state_full.getForces(asNumpy=True).value_in_unit(openmm.unit.kilojoules_per_mole / openmm.unit.nanometer)

        # Viriale scalare
        virial = -np.sum(positions * forces) / 3.0  # kJ/mol

        ke = state_full.getKineticEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)

        # Volume in nm³ → converti in m³ per i conti SI, poi torna in bar
        V_nm3 = Lx * Ly * Lz  # nm³
        V_m3 = V_nm3 * 1e-27  # nm³ → m³

        # N_A = 6.022e23 mol⁻¹
        # 1 kJ/mol / nm³ = 1e3 J / (6.022e23 * 1e-27 m³)
        #                 = 1e3 / (6.022e-4) Pa ≈ 1.661e6 Pa = 16.61 bar
        conversion = 1.66054e1  # kJ/(mol·nm³) → bar

        pressure = (2.0 * ke + 3.0 * virial) / (3.0 * V_nm3) * conversion

        values += [f'{Lx:.4f}', f'{Ly:.4f}', f'{Lz:.4f}', f'{pressure:.4f}']
        return values """

class BoxReporter:
    """Reporter che aggiunge Lx, Ly, Lz al CSV."""

    def __init__(self, file, reportInterval, barostat_index: int = None):
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        self._barostat_index = barostat_index
        self._file.write("\"Step\",Lx (nm),Ly (nm),Lz (nm),Pressure (bar)\n")

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        
        box = state.getPeriodicBoxVectors(asNumpy=True)
        Lx = box[0, 0].value_in_unit(openmm.unit.nanometer)
        Ly = box[1, 1].value_in_unit(openmm.unit.nanometer)
        Lz = box[2, 2].value_in_unit(openmm.unit.nanometer)
        step = simulation.currentStep

        pressure = float('nan')
        if self._barostat_index is not None:
            barostat = simulation.system.getForce(self._barostat_index)
            pressure_vector = barostat.computeCurrentPressure(simulation.context)
            pressure = (pressure_vector[0] + pressure_vector[1] + pressure_vector[2]).value_in_unit(openmm.unit.bar) / 3

        self._file.write(f"{step},{Lx:.4f},{Ly:.4f},{Lz:.4f},{pressure:.4f}\n")
        self._file.flush()

    def __del__(self):
        self._file.close()

#Funzione per trovare il percorso della cartella di lavoro
def get_base_dir() -> Path:
    # Script .py normale
    if '__file__' in dir():
        return Path(__file__).parent
    # Jupyter / IPython
    else:
        return Path.cwd()

def save_outputs(base_dir: Path, output_dir: str = None):
    #Sposta i file di output in una cartella dedicata

    if output_dir is None:
        #Cartella con timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = base_dir / f"run_{timestamp}"
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)

    files = [
        "topology.pdb",
        "trajectory.pdb",
        "data_full.csv",
        "data.csv",
        "box.csv",
        "final_state.xml",
        "interchange.json",
        "run.log"
    ]

    for filename in files:
        src = base_dir / filename
        if src.exists():
            shutil.move(str(src), str(output_dir / filename))
            logging.info(f"Spostato: {filename} > {output_dir}")
        else:
            logging.warning(f"File {filename} non trovato. Saltato.")
    
    src = base_dir / "parameters.toml"
    if src.exists():
        shutil.copy(str(src), str(output_dir / "parameters.toml"))
        logging.info(f'Copiato: parameters.toml > {output_dir}')
    else: 
        logging.warning(f"File parameters.toml non trovato. Saltato.")
    
    logging.info(f"Output salvati in {output_dir}")

    lastRun_dir = get_base_dir() / "Last_Run"
    if lastRun_dir.exists():
        shutil.rmtree(lastRun_dir)
    lastRun_dir.mkdir()

    shutil.copytree(output_dir, lastRun_dir, dirs_exist_ok = True)

    return output_dir

#Definizione di una funzione per il logging di tutto il programma
logging.basicConfig(
    level = logging.INFO,
    format = '%(asctime)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(str(get_base_dir() / "run.log"))
    ]
)

def create_simulation(
    interchange: Interchange,
    pdb_stride: int = 500,
    trajectory_name: str = str(get_base_dir() / "trajectory.pdb"),
    cutoff: float = 1.0, #nm -- AGGIUNTO 4/3/26
) -> openmm.app.Simulation:

    #AGGIUNTO 4/3/26
    #interchange = set_minimal_box(interchange, cutoff = cutoff*unit.nanometer)

    # Vecchio
    integrator = openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond, #Coupling del sistema al termostato (1/ps)
        1 * openmm.unit.femtoseconds, #Time stride. Base = 1 fs
    )

    # barostat = openmm.MonteCarloBarostat(
    #     1.0 * openmm.unit.bar,
    #     293.15 * openmm.unit.kelvin,
    #     25,
    # )

    # Barostato MonteCarlo anisotropico -- 5/3/26

    barostat = openmm.MonteCarloAnisotropicBarostat(
        (1.0, 1.0, 1.0) * openmm.unit.bar,
        293.15 * openmm.unit.kelvin,
        True,   #Spostamento lungo x
        True,   #Spostamento lungo y
        False,  #Spostamento lungo z -- BLOCCATO
        25,
    )

    simulation = interchange.to_openmm_simulation(
        combine_nonbonded_forces=True,
        integrator=integrator,
        additional_forces=[barostat],
    )

    # Vecchio
    simulation.minimizeEnergy()    

    simulation.context.setVelocitiesToTemperature(300 * openmm.unit.kelvin)
    simulation.context.computeVirtualSites()

    barostat_index = None
    for i in range(simulation.system.getNumForces()):
        if isinstance(simulation.system.getForce(i), openmm.MonteCarloAnisotropicBarostat):
            barostat_index = i
            break
    
    if barostat_index is None:
        logging.warning("Barostato non trovato nel sistema: la pressione non verrà calcolata.")

    pdb_reporter = openmm.app.PDBReporter(trajectory_name, pdb_stride)
    state_data_reporter = openmm.app.StateDataReporter(
        str(get_base_dir() / "data.csv"),
        10,
        step = True,
        potentialEnergy=True,
        kineticEnergy = True,
        totalEnergy = True,
        temperature = True,
        volume = True,
        density=True,
    )
    box_reporter = BoxReporter(
        str(get_base_dir() / "box.csv"),
        10,  # stesso intervallo del StateDataReporter
        barostat_index = barostat_index,
    )
    
    simulation.reporters.append(box_reporter)
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(state_data_reporter)

    # thermo = ReporterEsteso(
    #     str(get_base_dir() / "thermo.csv"),
    #     10,
    #     step = True,
    #     potentialEnergy = True,
    #     kineticEnergy = True,
    #     totalEnergy = True,
    #     temperature = True,
    #     volume = True,
    #     density = True,
    # )
    # simulation.reporters.append(thermo)

    return simulation

def run_simulation(
    simulation: openmm.app.Simulation,
    n_steps: int = 5000
):
    print("Starting simulation")
    logging.info("Starting simulation")
    start_time = time.process_time()

    print("Step, Dimensioni sistema (nm), volume (nm^3)")

    for step in range(n_steps):
        simulation.step(1)
        if step % 500 == 0:
            box_vectors = simulation.context.getState().getPeriodicBoxVectors()
            dimensioni_sistema = [box_vectors[i][i]._value for i in range(3)]
            print(step, dimensioni_sistema, np.linalg.det(box_vectors._value).round(3))


    end_time = time.process_time()
    print(f"Elapsed time: {(end_time - start_time):.2f} seconds")
    logging.info(f"Elapsed time: {(end_time - start_time):.2f} seconds")
    logging.info(f"Saving simulation")
    simulation.saveState(str(get_base_dir() / "final_state.xml"))
    logging.info("Simulation saved into final_state.xml")

def stima_box_volume(
    molecules: list, #lista molecole openFF
    n_copies: list, #nello stesso ordine della lista di molecole
    target_density: float, #g/mL
    z_fisso: float = None, #in nm. Se vuole bloccare Lz (voglio bloccarlo)
) -> np.ndarray:
#Stima dei vettori della box a partire da una densità target
#Meglio che z sia fisso.
#Restituisce il box in angstrom
    massa_totale = sum(
        sum(atom.mass.m for atom in mol.atoms) * n
        for mol, n in zip(molecules, n_copies)
    )

    N_a = 6.022e23 
    Volume_cella = (massa_totale / N_a) / (target_density * 1e-24)#in angstrom^3

    if z_fisso is not None:
        Lz = z_fisso * 10 # da nm a A
        Lxy = np.sqrt(Volume_cella / Lz)
        Lx = Lxy
        Ly = Lxy
    else:
        L = Volume_cella ** (1/3)
        Lx = Ly = Lz = L
    
    logging.info(f"Densità target: {target_density} g/mL")
    logging.info(f"Box: Lx = {Lx:.2f}, Ly = {Ly:.2f}, Lz = {Lz:.2f}")

    return np.array([
        [Lx, 0.0, 0.0],
        [0.0, Ly, 0.0],
        [0.0, 0.0, Lz]
    ]) * unit.angstrom

def build_from_toml(
    toml_parameters: dict,
    ):

    #Recupero molecole dalla directory fornita
    logging.info("Fetching molecules from given directories")

    HBD = Molecule.from_file(toml_parameters["molecules_paths"]["HBD_path"])
    HBA = Molecule.from_file(toml_parameters["molecules_paths"]["HBA_path"])
    controione = Molecule.from_file(toml_parameters["molecules_paths"]["counterion_path"])

    #Assegno residue names per tenere conto poi in visualizzazione o calcoli posteriori
    logging.info("Assigning residue names")

    for atom in HBD.atoms:
        atom.metadata["residue_name"] = "HBD"

    for atom in HBA.atoms:
        atom.metadata["residue_name"] = "HBA"

    for atom in controione.atoms:
        atom.metadata["residue_name"] = "CI"

    water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")

    for atom in water.atoms:
        atom.metadata["residue_name"] = "HOH"

    #Stessa cosa, ma per le due molecolone di cellobiosio
    logging.info("Fetching cellobiose unit cell")

    rdmol = Chem.MolFromPDBFile(toml_parameters["molecules_paths"]["cellobiose_topology"], removeHs = False)
    fragments = Chem.GetMolFrags(rdmol, asMols = True)

    #Separo le due catene
    offmol = [
        Molecule.from_rdkit(frag, allow_undefined_stereo = True)
        for frag in fragments
    ]

    #Assegno la carica formale alle due catene e assegno residue names diversi per le due catene (calcolerò la distanza?)

    i = 0
    for molecule in offmol:
        logging.info(("Assigning partial charge and residue names to chain:  " + str(i)))
        
        for atom in molecule.atoms:
            atom.metadata["residue_name"] = "CB"+str(i) #magari più avanti voglio misurare la distanza tra atomi dei due residui separati
        molecule.assign_partial_charges(toml_parameters["openff_methods"]["partial_charges_method"], use_conformers=molecule.conformers)
        i += 1
    polymer_top = Topology.from_molecules(offmol)

    polymer_positions = offmol[0].conformers[0]
    logging.info("Building solvent box")
    '''
    #dimensioni reali del polimero
    coords = polymer_positions.m_as("angstrom")

    xmin, ymin, zmin = coords.min(axis = 0)
    xmax, ymax, zmax = coords.max(axis = 0)

    Lx_poly = xmax - xmin
    Ly_poly = ymax - ymin
    Lz_poly = zmax - zmin

    #aggiungo padding sopra (x) e al lato (y)
    padding = toml_parameters["openff_methods"]["box_padding"] #angstrom
    Lx = (Lx_poly + padding)
    Ly = (Ly_poly + padding)
    Lz = Lz_poly

    # #definisco box
    box_vectors = np.array([
        [Lx, 0.0, 0.0],
        [0.0, Ly, 0.0],
        [0.0, 0.0, Lz]
    ]) * unit.angstrom
    '''
    polymer_pdb = openmm.app.PDBFile(toml_parameters["molecules_paths"]["cellobiose_topology"])
    cryst_box = polymer_pdb.topology.getPeriodicBoxVectors() #In nanometri

    #costruzione del box a partire da una densità stabilita (più robusto quando cambio le specie del DES)
    box_vectors = stima_box_volume(
        molecules = [HBD, HBA, controione, water],
        n_copies = [
            toml_parameters["molecules_copies"]["HBD_copies"],
            toml_parameters["molecules_copies"]["HBA_copies"],
            toml_parameters["molecules_copies"]["counterion_copies"],
            toml_parameters["molecules_copies"]["water_copies"],
        ],
        target_density = 1.2,
        z_fisso = cryst_box[2][2].value_in_unit(openmm.unit.nanometer),
    )
    # box_vectors = np.array(
    #     [
    #         [cryst_box[0][0].value_in_unit(openmm.unit.angstrom), 0.0, 0.0],
    #         [0.0, cryst_box[1][1].value_in_unit(openmm.unit.angstrom), 0.0],
    #         [0.0, 0.0, cryst_box[2][2].value_in_unit(openmm.unit.angstrom)]
    #         ]) * unit.angstrom
    print(type(box_vectors))
    print(box_vectors)

    # Solvente con pack_box. Include il polimero
    solvent_top = pack_box(
        molecules=[HBD, HBA, controione, water],
        number_of_copies=[
            toml_parameters["molecules_copies"]["HBD_copies"],
            toml_parameters["molecules_copies"]["HBA_copies"],
            toml_parameters["molecules_copies"]["counterion_copies"],
            toml_parameters["molecules_copies"]["water_copies"]
        ],
        box_vectors = box_vectors,
        solute = polymer_top
    )
    #Costruisco topologia definitiva
    logging.info("Building topology")

    final_topology = solvent_top
    final_topology.box_vectors = box_vectors

    final_topology.to_file(str(get_base_dir() / "topology.pdb"))

    #Costruisco forcefield e interchange openMM
    logging.info("Defining forcefield")

    sage = ForceField("openff-2.0.0.offxml")

    logging.info("BUilding interchange")
    interchange: Interchange = Interchange.from_smirnoff(force_field=sage, topology=final_topology)
    print("Interchange \nn_atoms:        " + str(interchange.topology.n_atoms) + "\nBox:            " + str(interchange.box) + "\nShape:          " + str(interchange.positions.shape))

    logging.info("Interchange \nn_atoms:        " + str(interchange.topology.n_atoms) + "\nBox:            " + str(interchange.box) + "\nShape:          " + str(interchange.positions.shape))
    
    #Salva interchange in un file appropriato
    open(str (get_base_dir() / "interchange.json"), 'w').write(interchange.model_dump_json())
    
    #Simulazione vera e propria
    logging.info("Creating simulation")
    simulation = create_simulation(interchange)

    return simulation

def build_from_state(
    toml_parameters: dict
):
    interchange = Interchange.model_validate_json(
        open(toml_parameters["resume_simulation"]["interchange_path"]).read()
    )
        
    simulation = create_simulation(interchange)
    simulation.loadState(toml_parameters["resume_simulation"]["state_path"])

    return simulation

def merge_csv(
    path1: str = str(get_base_dir() / "data.csv"),
    path2: str = str(get_base_dir() / "box.csv")
):
    #Prima ho creato due file csv con i due reporter separati. Adesso unisco i file csv in un unico file e cancello i vecchi
    logging.info("Creating CSV file")
    pandadata = pd.read_csv(path1)
    pandabox = pd.read_csv(path2)

    pandadata = pandadata.rename(columns={"#\"Step\"": "Step"})
    pandadata['Enthalpy (kJ/mole)'] = (pandadata["Total Energy (kJ/mole)"] * 1e3 + pandadata["Box Volume (nm^3)"] * 6.022e1 * pandabox["Pressure (bar)"]) * 1e-3
    merged = pd.merge(pandadata, pandabox, on = "Step")
    nuovo_ordine = [
        "Step",
        "Temperature (K)",
        "Potential Energy (kJ/mole)",
        "Kinetic Energy (kJ/mole)",
        "Total Energy (kJ/mole)",
        "Enthalpy (kJ/mole)",
        "Lx (nm)",
        "Ly (nm)",
        "Lz (nm)",
        "Pressure (bar)",
        "Box Volume (nm^3)",
        "Density (g/mL)"
    ]
    merged = merged[nuovo_ordine]
    merged.to_csv(str(get_base_dir() / "data_full.csv"))
    return 0


def main():
    #Programma per simulazioni di dinamica molecolare. Può prendere in input un file .toml di configurazione per eseguire un calcolo ab initio oppure un file openmm.Simulation.saveState per riprendere una simulazione precedente
    
    #Parser di argomenti.
    parser = argparse.ArgumentParser(
        prog = 'TestRun',
        description = 'Simulazioni di dinamica molecolare ab initio o da openMM.Simulation.loadState',
        add_help = True,
    )
    parser.add_argument('-f', type = str, required = True, default = str(get_base_dir() / "parameters.toml"), help="Percorso del file parameters.toml di configurazione della simulazione, sia che sia ab initio che una ripresa.")
    parser_group = parser.add_mutually_exclusive_group(required = True) #il calcolo ab initio e la ripresa sono mutualmente esclusivi
    parser_group.add_argument('--from_toml', type = bool, default = False)
    parser_group.add_argument('--from_state', type = bool, default = False)

    args = parser.parse_args()

    #Caricamento file TOML di configurazione
    logging.info("Reading parameters.toml file")
    with open(args.f, "rb") as f:
        toml_parameters = tomllib.load(f)
    
    with open(args.f, "rb") as f:
        toml_parameters = tomllib.load(f)
    
    #Casi specifici a seconda del parser degli argomenti forniti
    if args.from_toml:
        logging.info("Building simulation from scratch")
        simulation = build_from_toml(toml_parameters)
        steps = toml_parameters["openff_methods"]["simulation_steps"]

    elif args.from_state:
        logging.info("Resuming simulation")
        simulation = build_from_state(toml_parameters)
        steps = toml_parameters["resume_simulation"]["simulation_steps"]
    else:
        logging.info("NO ARGUMENT")
        return 1
    
    logging.info("Running simulation")
    run_simulation(simulation, steps)
    
    merge_csv()

    save_outputs(get_base_dir())
    logging.info("Done.")
    return 0

if __name__ == "__main__":
    main()
