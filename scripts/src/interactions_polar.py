import sys
from pathlib import Path
import mdtraj as md
import subprocess
import logging
import pickle
from collections import namedtuple, defaultdict
Chains = namedtuple('Chains', ['antibody', 'antigen'])
HBondAtom = namedtuple('HBondAtom', ['chainID', 'chain_type',
                       'CDR', 'resSeq', 'resname', 'index', 'serial', 'element', 'is_sidechain'])
HBond = namedtuple('HBond', ['donor', 'acceptor'])
# Water mediated doesnt' care about donor/acceptor.
WBond = namedtuple('WBond', ['water', 'residue'])

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)

AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")


def is_salt_bridge(donor, acceptor):
    donor_anion = (donor.resname == 'ASP' or donor.resname == 'GLU')\
        and donor.is_sidechain
    acceptor_anion = (acceptor.resname == 'ASP' or acceptor.resname == 'GLU')\
        and acceptor.is_sidechain
    donor_cation = (donor.resname == 'LYS' or donor.resname == 'ARG')\
        and donor.is_sidechain
    acceptor_cation = (acceptor.resname == 'LYS' or acceptor.resname == 'ARG')\
        and acceptor.is_sidechain

    return ((donor_anion and acceptor_cation) or (donor_cation and acceptor_anion))


def water_donor(linea, chains):

    donor_is_wat = linea[6:9] == 'HOH'
    resname = linea[20:23]
    chainID = linea[14]
    is_aminoacid = resname in AA_LIST
    is_in_ag = chainID in chains

    return (donor_is_wat and is_aminoacid and is_in_ag)


def water_acceptor(linea, chains):

    acceptor_is_wat = linea[20:23] == 'HOH'
    resname = linea[6:9]
    chainID = linea[0]
    is_aminoacid = resname in AA_LIST
    is_in_ag = chainID in chains

    return (acceptor_is_wat and is_aminoacid and is_in_ag)


def hbplus_lines_to_set(wat_lines, role_of_water="donor"):
    wat_dict = {}
    don_range = slice(0, 14)
    acc_range = slice(14, 28)

    if role_of_water == 'donor':
        wat_range = don_range
        mol_range = acc_range
    elif role_of_water == 'acceptor':
        wat_range = acc_range
        mol_range = don_range
    else:
        raise ValueError(f"hbplus_lines_to_set(): {role_of_water} is invalid argument "
                         f"to role_of_water")

    for linea in wat_lines:
        wat = linea[wat_range]
        if wat in wat_dict:
            wat_dict[wat].append(linea[mol_range])
        else:
            wat_dict[wat] = [linea[mol_range]]

    return wat_dict


def is_interface(linea, ab_chains, ag_chains):
    ag_donor_ab_acceptor = linea[0] in ag_chains and linea[14] in ab_chains
    ab_donor_ag_acceptor = linea[0] in ab_chains and linea[14] in ag_chains

    resname_donor = linea[6:9]
    resname_acceptor = linea[20:23]
    both_aminoacid = (resname_donor in AA_LIST) and (resname_acceptor in AA_LIST)

    return ((ag_donor_ab_acceptor or ab_donor_ag_acceptor) and both_aminoacid)


def get_hbonding_atom(topologia, chainID, resSeq, name):
    for cadena in topologia.chains:
        if cadena.chain_id == chainID:
            for resi in cadena.residues:
                if resi.resSeq == resSeq:
                    for atm in resi.atoms:
                        if atm.name == name:
                            return atm.index, atm.serial,\
                                atm.element.symbol, atm.is_sidechain
    raise ValueError('Hbonding atom was not found. Probably not an actual Hbond.')


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


def get_atom_from_string(pdb_idcode, linea, topologia, df_dataset, ab_chains):
    _chainID = linea[0]
    _resSeq = int(linea[1:5])
    _resname = linea[6:9]
    _chain_type, _cdr = get_chain_info(
        df_dataset, pdb_idcode, ab_chains, _chainID, _resSeq)
    _name = linea[10:14].strip()

    _index, _serial, _element, _is_sidechain = get_hbonding_atom(
        topologia, _chainID, _resSeq, _name)

    return HBondAtom(chainID=_chainID, resSeq=_resSeq, resname=_resname,
                     chain_type=_chain_type, CDR=_cdr, index=_index, serial=_serial,
                     element=_element, is_sidechain=_is_sidechain)


def get_interfacing_waters(wat_donor_ag_lines, wat_donor_ab_lines,
                           wat_acceptor_ag_lines, wat_acceptor_ab_lines):

    # Tried really hard to isolate donor/acceptor pairs, but waters can
    # act as double donor(1 to the antigen, another one to the antibody),
    # so it becomes even harder to keep track of this, and it turns up we
    # don't even care about this, so now `WBond` tuple is different from `Hbond`
    # tuple and there isn't a notion of donor/acceptor but one of water/residue
    wat_donor_ag_dict = hbplus_lines_to_set(wat_donor_ag_lines, "donor")
    wat_donor_ab_dict = hbplus_lines_to_set(wat_donor_ab_lines, "donor")
    wat_acceptor_ag_dict = hbplus_lines_to_set(wat_acceptor_ag_lines, "acceptor")
    wat_acceptor_ab_dict = hbplus_lines_to_set(wat_acceptor_ab_lines, "acceptor")

    # # Compile instances where water is the donor
    # wat_donor = wat_donor_ag_dict | wat_donor_ab_dict
    # interface_don_ag_acc_ab = set(wat_donor_ag_dict.keys()) &\
    #     set(wat_acceptor_ab_dict.keys())
    # interface_don_ag_don_ab = set(wat_donor_ag_dict.keys()) &\
    #     set(wat_donor_ab_dict.keys())
    # interface_don = {}
    # for wat in [*interface_don_ag_acc_ab, *interface_don_ag_don_ab]:
    #     interface_don[wat] = wat_donor[wat]

    # # Compile instances where water is the acceptor
    # wat_acceptor = wat_acceptor_ab_dict | wat_acceptor_ag_dict
    # interface_acc_ag_don_ab = set(wat_acceptor_ag_dict.keys())\
    #     & set(wat_donor_ab_dict.keys())
    # interface_acc_ag_acc_ab = set(wat_acceptor_ag_dict.keys())\
    #     & set(wat_acceptor_ab_dict.keys())
    # interface_acc = {}
    # for wat in [*interface_acc_ag_don_ab, *interface_acc_ag_acc_ab]:
    #     interface_acc[wat] = wat_acceptor[wat]

    wat_ag = wat_donor_ag_dict | wat_acceptor_ag_dict
    wat_ab = wat_donor_ab_dict | wat_acceptor_ab_dict
    interfacing_waters_set = set(wat_ag.keys()) & set(wat_ab.keys())

    interfacing_waters_ag = {wat: wat_ag[wat] for wat in interfacing_waters_set}
    interfacing_waters_ab = {wat: wat_ab[wat] for wat in interfacing_waters_set}

    return interfacing_waters_ag, interfacing_waters_ab


def parse_hb2(pdb_idcode, hb2_file, df_dataset, topologia, ab_chains, ag_chains):
    hbonds_dict = {}
    saltBridge_dict = {}
    wat_dict = defaultdict(list)

    wat_donor_ag_lines = []
    wat_acceptor_ag_lines = []
    wat_donor_ab_lines = []
    wat_acceptor_ab_lines = []
    with open(hb2_file, 'r') as archivo:

        lineas = archivo.readlines()
        for i, linea in enumerate(lineas):
            if i < 8:
                continue
            if water_donor(linea, ag_chains):
                wat_donor_ag_lines.append(linea)
            elif water_donor(linea, ab_chains):
                wat_donor_ab_lines.append(linea)
            elif water_acceptor(linea, ag_chains):
                wat_acceptor_ag_lines.append(linea)
            elif water_acceptor(linea, ab_chains):
                wat_acceptor_ab_lines.append(linea)
            elif is_interface(linea, ab_chains, ag_chains):
                try:
                    donor = get_atom_from_string(
                        pdb_idcode, linea[0:14], topologia, df_dataset, ab_chains)
                    acceptor = get_atom_from_string(
                        pdb_idcode, linea[14:28], topologia, df_dataset, ab_chains)

                    if is_salt_bridge(donor, acceptor):
                        saltBridge_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                        saltBridge_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                    else:
                        if donor.serial in hbonds_dict:
                            hbonds_dict[donor.serial].append(HBond(donor=donor, acceptor=acceptor))
                        else:
                            hbonds_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                        if acceptor.serial in hbonds_dict:
                            hbonds_dict[acceptor.serial].append(
                                HBond(donor=donor, acceptor=acceptor))
                        else:
                            hbonds_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                except ValueError:
                    logging.warning(f"Could not parse hb2 line: {i} {linea}")
                    continue

    # Now deal with the waters
    # interface_don, interface_acc = get_interfacing_waters(
    #     wat_donor_ag_lines, wat_donor_ab_lines, wat_acceptor_ag_lines, wat_acceptor_ab_lines)

     # Now deal with the waters
    interfacing_waters_ag, interfacing_waters_ab = get_interfacing_waters(
        wat_donor_ag_lines, wat_donor_ab_lines, wat_acceptor_ag_lines, wat_acceptor_ab_lines)

    for wat_str, res_list in [*interfacing_waters_ag.items(),
                              *interfacing_waters_ab.items()]:
        for res_str in res_list:
            wat = get_atom_from_string(
                pdb_idcode, wat_str, topologia, df_dataset, ab_chains)
            resi = get_atom_from_string(
                pdb_idcode, res_str, topologia, df_dataset, ab_chains)
            wat_dict[wat.serial].append(WBond(water=wat, residue=resi))

            # # Only add the protein residue as key
            # if acceptor.serial in wat_dict:
            #   wat_dict[acceptor.serial].append(WBond(water=wat, residue=resi))
            # else:
            #     wat_dict[acceptor.serial] = [WBond(water=donor, residue=acceptor)]

    # for wat_str, res_list in interface_don.items():
    #     for res_str in res_list:
    #         donor = get_atom_from_string(
    #             pdb_idcode, wat_str, topologia, df_dataset, ab_chains)
    #         acceptor = get_atom_from_string(
    #             pdb_idcode, res_str, topologia, df_dataset, ab_chains)

    #         # Only add the protein residue as key
    #         if acceptor.serial in wat_dict:
    #             wat_dict[acceptor.serial].append(HBond(donor=donor, acceptor=acceptor))
    #         else:
    #             wat_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]

    # for wat_str, res_list in interface_acc.items():
    #     for res_str in res_list:
    #         donor = get_atom_from_string(
    #             pdb_idcode, res_str, topologia, df_dataset, ab_chains)
    #         acceptor = get_atom_from_string(
    #             pdb_idcode, wat_str, topologia, df_dataset, ab_chains)
    #         # Only add the protein residue as key
    #         if donor.serial in wat_dict:
    #             wat_dict[donor.serial].append(HBond(donor=donor, acceptor=acceptor))
    #         else:
    #             wat_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]

    return hbonds_dict, saltBridge_dict, wat_dict


if __name__ == "__main__":
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
    print("Starting now.")

    hbonds_dict = {}
    saltBridge_dict = {}
    wat_dict = {}

    # check_pdb = '6edu'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path.joinpath(str_dir, pdb_idcode + '.pdb')
        try:
            trj_in = md.load(pdb_filename)
        except Exception as e:
            logging.error(f" Couldn't read {pdb_idcode}. Skipping.")
            continue
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        hbond_dir = Path.joinpath(data_dir, "hbonds")
        hbplus = Path.joinpath(source_location, 'hbplus')

        process = subprocess.run([hbplus, pdb_filename, "-A", "0", "0", "0", "-d", "3.9"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, cwd=hbond_dir)

        hb2_file = Path.joinpath(hbond_dir, pdb_filename.name[0:-3] + "hb2")

        hbonds_dict[pdb_idcode], saltBridge_dict[pdb_idcode], wat_dict[pdb_idcode] =\
            parse_hb2(pdb_idcode, hb2_file, epitope_buried_cleaned,
                      trj_in.topology, ab_chains, ag_chains)

    with open(Path.joinpath(casa_dir, 'data', 'polar_hbonds.pkl'), 'wb') as file:
        pickle.dump(hbonds_dict, file)
    with open(Path.joinpath(casa_dir, 'data', 'polar_saltBridge.pkl'), 'wb') as file:
        pickle.dump(saltBridge_dict, file)
    with open(Path.joinpath(casa_dir, 'data', 'polar_wat.pkl'), 'wb') as file:
        pickle.dump(wat_dict, file)
