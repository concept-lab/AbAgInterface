import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
AA_count = namedtuple('AA_count', ['heavy', 'light', 'antigen'])
HBondAtom = namedtuple('HBondAtom', ['chainID', 'chain_type',
                       'CDR', 'resSeq', 'resname', 'index', 'serial', 'element', 'is_sidechain'])
HBond = namedtuple('HBond', ['donor', 'acceptor'])
PolarCount = namedtuple('PolarCount', ['cdr_SC', 'cdr_BB', 'epi_SC', 'epi_BB'])
source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)
from abag_interactions_rings import *

casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
if __name__ == '__main__':

    logging.basicConfig(filename="log_" + Path(__file__).name, level=logging.INFO)
    logging.info("\n\n #### START #### \n")

    print("Reading epitope_buried_cleaned.pickle and other pickled data.")
    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]
    print("Starting now.")

    count = {}
    # check_pdb = '6ss6'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chain = chains[pdb_idcode].antigen

        heavy_cnt = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        light_cnt = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        antig_cnt = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        count_pdb = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        # In `ab_chains`, heavy chain is first, then the light one.
        # This comes from the abag naming convention in structures/exposed PDBs
        for chain in trj_in.topology.chains:
            if chain.chain_id == ab_chains[0]:
                for resi in chain.residues:
                    heavy_cnt[resi.name[0:3]] += 1
            elif chain.chain_id == ab_chains[1]:
                for resi in chain.residues:
                    light_cnt[resi.name[0:3]] += 1
            elif chain.chain_id in ag_chain:
                for resi in chain.residues:
                    try:
                        antig_cnt[resi.name[0:3]] += 1
                    except KeyError:
                        logging.warning(f"Weird resname: {resi.name[0:3]} in pdb: {pdb_idcode}")

        for aa in count_pdb.keys():
            try:
                count_pdb[aa] = AA_count(
                    heavy=heavy_cnt[aa], light=light_cnt[aa], antigen=antig_cnt[aa])
            except KeyError:
                logging.warning(f"Weird resname: {aa[0:3]} in pdb: {pdb_idcode}")
        
        count[pdb_idcode] = count_pdb
        
    with open(Path.joinpath(casa_dir, 'data', 'count_aa.pkl'), 'wb') as file:
        pickle.dump(count, file)


  