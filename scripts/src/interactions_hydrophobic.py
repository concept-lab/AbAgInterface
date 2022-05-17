import sys
from pathlib import Path
import pickle
import pandas as pd
import mdtraj as md
import networkx as nx
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
InterfaceAtoms = namedtuple('InterfaceAtoms', ['antibody', 'antigen'])
Atom = namedtuple('Atom', ['index', 'serial', 'element', 'is_sidechain', 'resSeq',
                  'resSeq_str', 'resname', 'chain_ID', 'chain_type', 'CDR'])

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)
from abag_interactions_hydrophobic import *

casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

if __name__ == '__main__':
    logging.basicConfig(filename="log_" + Path(__file__).name, level=logging.INFO)

    print("Reading all the .pickle's")
    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]
    with open(Path.joinpath(data_dir, 'interface_atoms.pkl'), 'rb') as file:
        interface_atoms = pickle.load(file)
    print("Starting.")

    hydrophobic_clusters = {}
    check_pdb = '7jmp'
    idx = pdb_list.index(check_pdb)
    for pdb_idcode in [pdb_list[idx]]:
        # for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        try:
            G = get_carbons_graph(trj_in, interface_atoms[pdb_idcode], cutoff_carbons)
            pre_clusteres = get_putative_clusters(G)
            clusters = merge_clusters(trj_in, pre_clusteres, cutoff_clusters)
        except Exception as e:
            logging.error(
                f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during hydrophobic "
                f"interactions calculation. Probably has no hydrophobic interactions. "
                f"Assigning an empty tuple.")
            hydrophobic_clusters[pdb_idcode] = ()
            continue
        else:
            hydrophobic_clusters[pdb_idcode] = clusters

    # with open(Path.joinpath(data_dir, 'hydrophobic.pkl'), 'wb') as file:
    #     pickle.dump(hydrophobic_clusters, file)

    print(f" --- Done -- ")
