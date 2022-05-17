import sys
from pathlib import Path
import itertools
import pickle
import mdtraj as md
from collections import Counter
from collections import namedtuple
ResiCount = namedtuple('ResiCount', ['antibody', 'antigen'])
PolarCount = namedtuple('PolarCount', ['cdr_SC', 'cdr_BB', 'epi_SC', 'epi_BB'])
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
    with open(Path.joinpath(casa_dir, "data", 'epitope_buried_cleaned.pickle'), 'rb') as file:
        epitope_buried_cleaned = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]
    print("Starting now.")
    
    count_resi_inteface_pdb = {}
    count_resi_inteface = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
    count_tot_pdb = {}
    count_tot_chain_pdb = {}
    count_tot = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)
        
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))

        count_res_ab = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        count_res_ag = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))

        ###
        # TOTAL
        ###
        count = Counter([residue.name for residue in trj_in.topology.residues])
        count_tot_this_pdb = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        for key, val in count.items():
            if key not in AA_LIST:
                continue
            count_tot_this_pdb[key] = val
            count_tot[key] += val
        count_tot_pdb[pdb_idcode] = count_tot_this_pdb

        cnt_ag = len([residue for residue in trj_in.topology.residues\
            if residue.chain.chain_id in chains[pdb_idcode].antigen ])
        cnt_ab = len([residue for residue in trj_in.topology.residues\
            if residue.chain.chain_id in chains[pdb_idcode].antibody ])
        count_tot_chain_pdb[pdb_idcode] = ResiCount(antibody=cnt_ab, antigen=cnt_ag)
        ###
        # ANTIBODY INTERFACE
        ###
        interface_ab = []
        for linea in epitope_buried_cleaned.query(f"idcode == '{pdb_idcode}'").ab_ag_interface_res.values[0]:
            chainID = linea[0]
            resSeq_str = linea[1]
            resname = linea[2]
            # Including resname in the mix due to heavychain CDR3 extra residues with
            # a letter at the end. They will show up with the same resSeq in mdtraj
            # topology, but as long as I'm querying them, there's no problem.
            interface_ab.append(chainID+resname+resSeq_str)

        for unique_ in set(interface_ab):
            resname = unique_[1:4]
            count_res_ab[resname] += 1

        ###
        # ANTIGEN INTERFACE
        ###
        interface_ag = []
        for linea in epitope_buried_cleaned.query(f"idcode == '{pdb_idcode}'").ag_ab_interface_res.values[0]:
            chainID = linea[0]
            resSeq_str = linea[1]
            resname = linea[2]
            # Including resname in the mix due to heavychain CDR3 extra residues with
            # a letter at the end. They will show up with the same resSeq in mdtraj
            # topology, but as long as I'm querying them, there's no problem.
            interface_ag.append(chainID+resname+resSeq_str)

        for unique_ in set(interface_ag):
            resname = unique_[1:4]
            count_res_ag[resname] += 1

        count_resi_inteface_pdb[pdb_idcode] = ResiCount(antibody=count_res_ab, antigen=count_res_ag)
        for key in count_resi_inteface.keys():
            count_resi_inteface[key] += count_res_ab[key]
            count_resi_inteface[key] += count_res_ag[key]
        
    with open(Path.joinpath(casa_dir, 'data', 'count_resi_inteface_pdb.pkl'), 'wb') as file:
        pickle.dump(count_resi_inteface_pdb, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_resi_inteface.pkl'), 'wb') as file:
        pickle.dump(count_resi_inteface, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_tot_pdb.pkl'), 'wb') as file:
        pickle.dump(count_tot_pdb, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_tot_chain_pdb.pkl'), 'wb') as file:
        pickle.dump(count_tot_chain_pdb, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_tot.pkl'), 'wb') as file:
        pickle.dump(count_tot, file)



  