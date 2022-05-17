import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
InterfaceAtoms = namedtuple('InterfaceAtoms', ['antibody', 'antigen'])

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
    logging.info("\n\n #### START #### \n")

    print("Reading all the .pickle's")
    # with open(Path.joinpath(casa_dir, "data", 'epitope_buried_cleaned.pickle'), 'rb') as file:
    #     epitope_buried_cleaned = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]
    with open(Path.joinpath(data_dir, 'interface_atoms.pkl'), 'rb') as file:
        interface_atoms = pickle.load(file)
    print("Starting.")

    SSE = {}
    SSE_cnt = {}
    SSE_cnt_interface = {}
    # check_pdb = '1qfw'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ag_chains = chains[pdb_idcode].antigen

        # Get the resSeqs of the interface
        interface_resSeq = set(
            [atom.resSeq for atom in interface_atoms[pdb_idcode].antigen.values()])
        # Get SSE of each residue
        dssp = md.compute_dssp(trj_in, simplified=True)[0]
        SSE_pdb = {}
        SSE_cnt_pdb = {'H': 0, 'E': 0, 'C': 0, 'NA': 0}
        SSE_cnt_interface_pdb = {'H': 0, 'E': 0, 'C': 0, 'NA': 0}
        for i, res in enumerate(trj_in.topology.residues):
            if res.chain.chain_id in ag_chains:
                SSE_pdb[res.resSeq] = dssp[i]
                SSE_cnt_pdb[dssp[i]] += 1
                if res.resSeq in interface_resSeq:
                    SSE_cnt_interface_pdb[dssp[i]] += 1

        SSE[pdb_idcode] = SSE_pdb
        SSE_cnt[pdb_idcode] = SSE_cnt_pdb
        SSE_cnt_interface[pdb_idcode] = SSE_cnt_interface_pdb

    with open(Path.joinpath(casa_dir, 'data', 'simple_SSE.pkl'), 'wb') as file:
        pickle.dump(SSE, file)
    with open(Path.joinpath(casa_dir, 'data', 'simple_SSE_count.pkl'), 'wb') as file:
        pickle.dump(SSE_cnt, file)
    with open(Path.joinpath(casa_dir, 'data', 'simple_SSE_cnt_interface.pkl'), 'wb') as file:
        pickle.dump(SSE_cnt_interface, file)
