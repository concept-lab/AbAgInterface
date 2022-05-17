import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
PionPair = namedtuple('PionPair', ['ring', 'ion'])
Ring = namedtuple('Ring', ['indices', 'serials', 'resSeq', 'resSeq_str',
                           'resname', 'chain_ID', 'chain_type', 'CDR'])
Atom = namedtuple('Atom', ['index', 'serial', 'element', 'is_sidechain', 'resSeq',
                  'resSeq_str', 'resname', 'chain_ID', 'chain_type', 'CDR'])
InterfaceAtoms = namedtuple('InterfaceAtoms', ['antibody', 'antigen'])

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)
from abag_interactions_rings import *

casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

if __name__ == '__main__':
    logging.basicConfig(filename="log_" + Path(__file__).name, level=logging.INFO)
    logging.info("\n\n #### START #### \n")

    print("Reading epitope_buried_cleaned.pickle and other pickled data.")
    with open(Path.joinpath(
            casa_dir, "data", 'epitope_buried_cleaned.pickle'), 'rb') as file:
        epitope_buried_cleaned = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]
    with open(Path.joinpath(data_dir, 'interface_atoms.pkl'), 'rb') as file:
        interface_atoms = pickle.load(file)
    print("Starting now.")

    PiAnion = {}
    PiCation = {}
    # check_pdb = '5i5k'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen
        PiAnion[pdb_idcode] = tuple()
        PiCation[pdb_idcode] = tuple()

        try:
            CG_rings, CoM_rings_xyz, normal_vectors = get_ring_data(
                trj_in, ab_chains, ring_atoms_pi_ion)

            ab_anions, ag_anions, ab_cations, ag_cations = get_ions(
                interface_atoms[pdb_idcode])

            ring_ab_anion_ag = get_ion_ring_interactions(
                trj_in, CG_rings["antibody"], ag_anions, CoM_rings_xyz["antibody"],
                normal_vectors["antibody"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)
            ring_ab_cation_ag = get_ion_ring_interactions(
                trj_in, CG_rings["antibody"], ag_cations, CoM_rings_xyz["antibody"],
                normal_vectors["antibody"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)

            ring_ag_anion_ab = get_ion_ring_interactions(
                trj_in, CG_rings["antigen"], ab_anions, CoM_rings_xyz["antigen"],
                normal_vectors["antigen"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)

            ring_ag_cation_ab = get_ion_ring_interactions(
                trj_in, CG_rings["antigen"], ab_cations, CoM_rings_xyz["antigen"],
                normal_vectors["antigen"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)

            data_ring_ab_anion_ag = get_data_from_ring_ab_ion_ag(
                pdb_idcode, epitope_buried_cleaned, ring_ab_anion_ag)

            data_ring_ab_cation_ag = get_data_from_ring_ab_ion_ag(
                pdb_idcode, epitope_buried_cleaned, ring_ab_cation_ag)

            data_ring_ag_anion_ab = get_data_from_ring_ag_ion_ab(
                pdb_idcode, epitope_buried_cleaned, ring_ag_anion_ab)

            data_ring_ag_cation_ab = get_data_from_ring_ag_ion_ab(
                pdb_idcode, epitope_buried_cleaned, ring_ag_cation_ab)

        except Exception as e:
            logging.error(
                f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during Pi-ion "
                f"interactions calculation. Aborting.")
            raise e
        else:
            PiAnion[pdb_idcode] = (*data_ring_ab_anion_ag, *data_ring_ag_anion_ab)
            PiCation[pdb_idcode] = (*data_ring_ab_cation_ag, *data_ring_ag_cation_ab)

    with open(Path.joinpath(casa_dir, 'data', 'PiAnion.pkl'), 'wb') as file:
        pickle.dump(PiAnion, file)
    with open(Path.joinpath(casa_dir, 'data', 'PiCation.pkl'), 'wb') as file:
        pickle.dump(PiCation, file)

    print(f" --- Done -- ")
