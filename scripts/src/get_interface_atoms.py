import pickle
import sys
import string
from pathlib import Path
import logging
import mdtraj as md
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
InterfaceAtoms = namedtuple('InterfaceAtoms', ['antibody', 'antigen'])
Atom = namedtuple('Atom', ['index', 'serial', 'element', 'is_sidechain', 'resSeq',
                  'resSeq_str', 'resname', 'chain_ID', 'chain_type', 'CDR'])

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)

casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

# Get a reference to each of the mdtraj chains


def get_chains(topologia, ab_chains, ag_chains):
    cadenas = {}
    try:
        cadenas[ab_chains[0]] = next((chain for chain in topologia.chains
                                      if chain.chain_id == ab_chains[0]))
        cadenas[ab_chains[1]] = next((chain for chain in topologia.chains
                                      if chain.chain_id == ab_chains[1]))
        for chain_ in ag_chains:
            if chain_ != '.':
                # Only ~10% of the PDBs have more than 1 antigen chain.
                cadenas[chain_] = next((chain for chain in topologia.chains
                                        if chain.chain_id == chain_))
    except StopIteration:
        logging.error(f"{pdb_idcode} lacks one of these chains: {(ab_chains, ag_chains)}"
                      f"This shouldn't happen. Filename is: {filenames[pdb_idcode]}")
        raise
    return cadenas


def get_chain_types(epitope_buried_cleaned, pdb_idcode, ab_chains, ag_chains):
    chain_types = {}

    chain_types[ab_chains[0]] = epitope_buried_cleaned.query(
        f"idcode == '{pdb_idcode}' and chainID == '{ab_chains[0]}'").chain_type.values[0]
    assert chain_types[ab_chains[0]].isupper(), f"Invalid chain type: "
    f"{chain_types[ab_chains[0]]} for chainID: {ab_chains[0]}"

    chain_types[ab_chains[1]] = epitope_buried_cleaned.query(
        f"idcode == '{pdb_idcode}' and chainID == '{ab_chains[1]}'").chain_type.values[0]
    assert chain_types[ab_chains[1]].isupper(), f"Invalid chain type: "
    f"{chain_types[ab_chains[1]]} for chainID: {ab_chains[1]}"

    for chain_ in ag_chains:
        chain_types[chain_] = '.'

    return chain_types


def get_cdr_from_residue(epitope_buried_cleaned, pdb_idcode, chain_type, residue):
    # Antigen residue:
    cdr = -1
    if chain_type != '.':
        # Antibody residue:
        rows_for_this_chain = epitope_buried_cleaned.query(
            f"idcode == '{pdb_idcode}' and chainID == '{residue.chain.chain_id}'")
        for row in rows_for_this_chain.iterrows():
            try:
                cdr_begin = int(row[1].cdr_begin.strip())
            except ValueError:
                cdr_begin = int(row[1].cdr_begin.strip()[:-1])
            try:
                cdr_end = int(row[1].cdr_end.strip())
            except ValueError:
                cdr_end = int(row[1].cdr_end.strip()[:-1])

            if cdr_begin <= residue.resSeq <= cdr_end:
                cdr = row[1].cdr
                break
        else:
            # Framework residue:
            cdr = 0

    return cdr


def get_atoms_from_rows(pdb_idcode, epitope_buried_cleaned, atomic_lines, cadenas, chain_types):
    alfabeto = tuple([""] + list(string.ascii_uppercase))
    set_of_residues = set()
    atoms = {}

    for atom_line in atomic_lines:
        chainID = atom_line[0]
        resSeq_str = atom_line[1].strip()
        resname = atom_line[2]

        try:
            resSeq = int(resSeq_str)
        except ValueError:
            # Residue wth a character.
            # Once a residue with a character in its resSeq shows up, I deal
            # with all the following residues at once.
            resSeq = int(resSeq_str[:-1])
            unique_residue_id = chainID + resname + resSeq_str[:-1] + 'Z'
            if unique_residue_id in set_of_residues:
                continue
            set_of_residues.add(unique_residue_id)

            for i, residue in enumerate(
                    resi for resi in cadenas[chainID].residues if resi.resSeq == resSeq):
                if i == 0:
                    # The first residue is the one without a letter. Don't add it.
                    continue
                cdr = get_cdr_from_residue(
                    epitope_buried_cleaned, pdb_idcode, chain_types[chainID], residue)

                for atom in residue.atoms:
                    # I have to type name[0:3] because there was an actual fucker
                    # who used 4 characters for a residue name. Fuck that guy.
                    atoms[atom.serial] = Atom(
                        index=atom.index, serial=atom.serial, element=atom.element.symbol,
                        is_sidechain=atom.is_sidechain, resSeq=residue.resSeq,
                        resSeq_str=str(residue.resSeq) + alfabeto[i],
                        resname=residue.name[0:3], chain_ID=chainID,
                        chain_type=chain_types[chainID], CDR=cdr)

        else:
            # Regular residue.
            unique_residue_id = chainID + resname + resSeq_str
            if unique_residue_id in set_of_residues:
                continue
            set_of_residues.add(unique_residue_id)
            residue = next(resi for resi in cadenas[chainID].residues if resi.resSeq == resSeq)
            cdr = get_cdr_from_residue(
                epitope_buried_cleaned, pdb_idcode, chain_types[chainID], residue)
            for atom in residue.atoms:
                atoms[atom.serial] = Atom(
                    index=atom.index, serial=atom.serial, element=atom.element.symbol,
                    is_sidechain=atom.is_sidechain, resSeq=residue.resSeq, resSeq_str=resSeq_str,
                    resname=residue.name[0:3], chain_ID=chainID,
                    chain_type=chain_types[chainID], CDR=cdr)
    return atoms


if __name__ == '__main__':
    logging.basicConfig(filename="log_" + Path(__file__).name, level=logging.INFO)
    print("Reading buried_interface_res.pickle")
    with open(Path.joinpath(
            casa_dir, "data", 'buried_interface_res.pickle'), 'rb') as file:
        buried_interface_res = pickle.load(file)
    print("Reading epitope_buried_cleaned.pickle")
    with open(Path.joinpath(
            casa_dir, "data", 'epitope_buried_cleaned.pickle'), 'rb') as file:
        epitope_buried_cleaned = pickle.load(file)
    print("Starting now.")

    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'r')) as file:
        pdb_list = [linea.strip() for linea in file]

    interface_atoms = {}
    # check_pdb = '1ahw'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        cadenas = get_chains(trj_in.topology, ab_chains, ag_chains)
        chain_types = get_chain_types(epitope_buried_cleaned, pdb_idcode, ab_chains, ag_chains)

        # `buried_interface_res` has 1 row per interface. In case there's more than one,
        # I use the heavy chain (ab_chains[0]) to identify the one I'm care about.
        df_interface_atoms = buried_interface_res.query(f"idcode == '{pdb_idcode}'\
            and chainID == '{ab_chains[0]}'")

        try:
            ab_atoms = get_atoms_from_rows(
                pdb_idcode, epitope_buried_cleaned, df_interface_atoms.ab_ag_interface_res.values
                [0],
                cadenas, chain_types)
        except Exception as e:
            logging.exception(e)
            logging.error(f" {pdb_idcode} antibody's interface atom look-up failed.")
            continue

        try:
            ag_atoms = get_atoms_from_rows(
                pdb_idcode, epitope_buried_cleaned, df_interface_atoms.ag_ab_interface_res.values
                [0],
                cadenas, chain_types)
        except Exception as e:
            logging.exception(e)
            logging.error(f" {pdb_idcode} antigen's interface atom look-up failed.")
            continue

        interface_atoms[pdb_idcode] = InterfaceAtoms(antibody=ab_atoms, antigen=ag_atoms)

    with open(Path.joinpath(data_dir, 'interface_atoms.pkl'), 'wb') as file:
        pickle.dump(interface_atoms, file)
