import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
ShieldingAtom = namedtuple(
    'ShieldingAtom', ['chainID', 'chain_type', 'CDR', 'resSeq', 'resname', 'index',
                      'serial', 'element', 'is_sidechain'])
InterfaceAtoms = namedtuple('InterfaceAtoms', ['antibody', 'antigen'])
source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)
from abag_interactions_hydrophobic import *

AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")
cutoff = .5

#
# Using mdtraj IDs for shielding indices in the output dict.
# Should;ve used serials
#

def get_chain_info(epitope_buried_cleaned, pdb_idcode, ab_chains, chainID, resSeq):
    # Framework residue:
    cdr = 0
    if chainID in ab_chains:
        rows_for_this_chain = epitope_buried_cleaned.query(
            f"idcode == '{pdb_idcode}' and chainID == '{chainID}'")
        for row in rows_for_this_chain.iterrows():
            try:
                cdr_begin = int(row[1].cdr_begin.strip())
            except ValueError:
                cdr_begin = int(row[1].cdr_begin.strip()[:-1])
            try:
                cdr_end = int(row[1].cdr_end.strip())
            except ValueError:
                cdr_end = int(row[1].cdr_end.strip()[:-1])

            if cdr_begin <= resSeq <= cdr_end:
                cdr = row[1].cdr
                break

        chain_type = 'L' if rows_for_this_chain.chain_type.values[0] in (
            'K', 'L') else rows_for_this_chain.chain_type.values[0]
        assert chain_type in ('K', 'L', 'H'), f"{pdb_idcode} returned invalid "\
            f"chain type for hbonding residue: {resSeq=}, {chainID=}"
    else:
        chain_type = -1
        cdr = -1
    return chain_type, cdr


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

    shielding_dict = {}
    # check_pdb = '1adq'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        # pdb_filename = Path.joinpath(str_dir, pdb_idcode + '.pdb')
        pdb_filename = Path(filenames[pdb_idcode])
        try:
            trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
            # trj_in = md.load(pdb_filename)
        except Exception as e:
            logging.error(f" Couldn't read {pdb_idcode}. Skipping.")
            continue
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        ab_carbons = [atom.index for atom in interface_atoms[pdb_idcode].antibody.values()
                      if atom.element == 'C']

        ag_carbons = [atom.index for atom in interface_atoms[pdb_idcode].antigen.values()
                      if atom.element == 'C']

        ab_polars = [atom.index for atom in interface_atoms[pdb_idcode].antibody.values()
                     if atom.element in ('N', 'O')]
        ag_polars = [atom.index for atom in interface_atoms[pdb_idcode].antigen.values()
                     if atom.element in ('N', 'O')]
        polars = ab_polars + ag_polars

        C_ON_pairs = np.array(list(itertools.product(ab_carbons, polars)))
        C_C_pairs = np.array(list(itertools.product(ab_carbons, ag_carbons)))

        C_ON_distancias = md.compute_distances(
            trj_in, C_ON_pairs).reshape((len(ab_carbons), len(polars)))

        C_C_distancias = md.compute_distances(
            trj_in, C_C_pairs).reshape((len(ab_carbons), len(ag_carbons)))

        G = nx.Graph()
        indices_close_C_C_distancias = np.where(C_C_distancias < cutoff)
        mask_close_C_ON_distancias = C_ON_distancias < cutoff
        shielding_pdb = {}
        for i, j in zip(*indices_close_C_C_distancias):
            C_cdr_id = ab_carbons[i]
            C_epi_id = ag_carbons[j]
            surrounding_ON_ids = [
                polars[i]
                for i in np.where(mask_close_C_ON_distancias[i, :])[0]]

            shielded, ON_id = is_shielded(
                trj_in.xyz[0], C_cdr_id, C_epi_id, surrounding_ON_ids)

            if shielded:
                # I could rewrite this to get the data from `interface_atoms[pdb_idcode]`
                # but I'm in a bit of a hurry.
                chainID = trj_in.topology.atom(ON_id).residue.chain.chain_id
                resSeq = trj_in.topology.atom(ON_id).residue.resSeq
                resname = trj_in.topology.atom(ON_id).residue.name
                chain_type, cdr = get_chain_info(
                    epitope_buried_cleaned, pdb_idcode, ab_chains, chainID, resSeq)
                serial = trj_in.topology.atom(ON_id).serial
                element = trj_in.topology.atom(ON_id).element.symbol
                is_sidechain = trj_in.topology.atom(ON_id).is_sidechain

                # Compile all the info on this shielding polar atom.
                shielding_atom = ShieldingAtom(
                    chainID=chainID, chain_type=chain_type, CDR=cdr,
                    resSeq=resSeq, resname=resname, index=ON_id,
                    serial=serial, element=element, is_sidechain=is_sidechain)
                
                shielding_pdb[ON_id] = shielding_atom

        shielding_dict[pdb_idcode] = shielding_pdb

    with open(Path.joinpath(casa_dir, 'data', 'shielding.pkl'), 'wb') as file:
        pickle.dump(shielding_dict, file)

    print(f" --- Done -- ")
