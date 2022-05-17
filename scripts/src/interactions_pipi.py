import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
PiPiPair = namedtuple('PiPiPair', ['antibody', 'antigen'])
Ring = namedtuple('Ring', ['indices', 'serials', 'resSeq', 'resSeq_str',
                           'resname', 'chain_ID', 'chain_type', 'CDR'])

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
    print("Starting now.")

    PiPi = {}
    # check_pdb = '4oqt'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen
        PiPi[pdb_idcode] = tuple()

        try:
            CG_rings, CoM_rings_xyz, normal_vectors = get_ring_data(
                trj_in, ab_chains, ring_atoms)
            pipi_ring_pairs = get_pipi_interactions(
                trj_in, CG_rings, CoM_rings_xyz, normal_vectors, cutoff_ring,
                cutoff_angle_pipi)

            rings = get_data_from_ring_ring(pdb_idcode, epitope_buried_cleaned, pipi_ring_pairs)
        except Exception as e:
            logging.error(
                f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during PiPi "
                f"interactions calculation.")
            raise e
        else:
            PiPi[pdb_idcode] = rings

    with open(Path.joinpath(casa_dir, 'data', 'PiPi.pkl'), 'wb') as file:
        pickle.dump(PiPi, file)

    print(f" --- Done -- ")
