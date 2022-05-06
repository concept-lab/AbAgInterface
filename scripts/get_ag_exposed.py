import os
import pandas as pd
from utils import read_pdb_line, get_sabdab_details, pqr2xyzr

from get_exposed import iterate_prots, System, get_interfaces
from constants import EXPOSED_DIR, RAW_STRUCTURES_DIR


def get_systems(idcode, antigen_chains, cdr_chain, ab_chains):
    chain_res = {}
    # print(f"{RAW_STRUCTURES_DIR}/{idcode}.pdb")
    with open(f"{RAW_STRUCTURES_DIR}/{idcode}.pdb") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            (
                chain,
                cdr_id,
                (resname, resnumb),
                (aname, anumb),
                (x, y, z),
            ) = read_pdb_line(line)

            if chain not in chain_res:
                chain_res[chain] = []
            chain_res[chain].append(resnumb)

    if ab_chains[0] != ab_chains[0] or ab_chains[1] != ab_chains[1]:
        otherab_chain = ""
        chains = ab_chains[0] if ab_chains[1] != ab_chains[1] else ab_chains[1]
    else:
        chains = "".join(ab_chains)
        otherab_chain = chains.replace(cdr_chain, "").strip()

    f_complex = f"{idcode}_complex_{chains}_{antigen_chains}"
    f_antibody = f"{idcode}_antibody_{chains}"
    f_antigen = f"{idcode}_antigen_{antigen_chains}"

    prot_dir = f"{EXPOSED_DIR}/{idcode}/"

    ab_ag_complex = System(
        f_complex,
        {
            chain: list(set(chain_res[chain]))
            for chain in antigen_chains + cdr_chain + otherab_chain
        },
        prot_dir,
    )
    antibody = System(
        f_antibody,
        {chain: list(set(chain_res[chain])) for chain in cdr_chain + otherab_chain},
        prot_dir,
    )
    antigen = System(
        f_antigen,
        {chain: list(set(chain_res[chain])) for chain in antigen_chains},
        prot_dir,
    )

    return (ab_ag_complex, antibody, antigen)


if __name__ == "__main__":
    df_sabdab_all = get_sabdab_details()

    f_epitopes = "../data/epitope_buried.pickle"  # "../data/cdr_epitope.pickle"

    ab_int_aa_dist = {}
    ag_int_aa_dist = {}
    ag_aa_dist = {}
    for prot in iterate_prots(f_epitopes):
        if prot["idcode"] == "6kn9":
            continue

        if not prot["ab_ag_interface_res"]:
            continue

        ab_chain_type = prot["chain_type"]
        if prot["chain_type"] == "K":
            ab_chain_type = "L"
        df_sabdab_prot = df_sabdab_all.query(
            f"pdb == '{prot['idcode']}' and {ab_chain_type}chain == '{prot['chainID']}'"
        )

        antigen_chains = df_sabdab_prot["antigen_chain"].values[0].replace(" | ", "")
        antibody_chains = df_sabdab_prot[["Hchain", "Lchain"]].values[0]

        prot_dir = f'{EXPOSED_DIR}/{prot["idcode"]}/'

        (ab_ag_complex, antibody, antigen) = get_systems(
            prot["idcode"],
            antigen_chains,
            prot["chainID"],
            antibody_chains,
        )

        atom_res = pqr2xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr, prot["cdr"])

        for system in (ab_ag_complex, antibody, antigen):
            system.trim_xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr)
            if os.path.isfile(system.fexposed):
                system.get_surface_atoms()

        ag_ab_interface = antigen.surface_atoms - ab_ag_complex.surface_atoms
        ag_ab_interface_res = [atom_res[i][:3] for i in ag_ab_interface]

        ab_ag_interface = antibody.surface_atoms - ab_ag_complex.surface_atoms
        ab_ag_interface_res = [atom_res[i][:3] for i in ab_ag_interface]

        ag_interface = antigen.surface_atoms
        ag_interface_res = [atom_res[i][:3] for i in ag_interface]

        dists = [
            (ab_ag_interface_res, ab_int_aa_dist),
            (ag_ab_interface_res, ag_int_aa_dist),
            (ag_interface_res, ag_aa_dist),
        ]
        for residues, dist in dists:
            for res in set(residues):
                aa = res[2]
                if aa not in dist:
                    dist[aa] = 0
                dist[aa] += 1

    print(ab_int_aa_dist)
    print(ag_int_aa_dist)
    print(ag_aa_dist)


# {'ASN': 27660, 'THR': 28146, 'SER': 55002, 'TRP': 26361, 'GLY': 34839, 'ARG': 23160, 'PHE': 19935, 'TYR': 84096, 'ILE': 14724, 'ASP': 23655, 'GLU': 13578, 'HIS': 8793, 'VAL': 15270, 'LYS': 8277, 'LEU': 17661, 'GLN': 7758, 'MET': 2664, 'ALA': 17238, 'PRO': 11190, 'CYS': 1248}
# {'GLU': 24051, 'PRO': 22338, 'ARG': 25332, 'LEU': 28845, 'HIS': 10167, 'ILE': 20958, 'SER': 25839, 'PHE': 15579, 'GLN': 20178, 'THR': 31017, 'TRP': 8049, 'CYS': 6369, 'ALA': 19938, 'MET': 5652, 'VAL': 20280, 'TYR': 22338, 'GLY': 31086, 'ASN': 30438, 'LYS': 28206, 'ASP': 25218}
# {'GLU': 432135, 'ALA': 447894, 'LEU': 622368, 'THR': 582390, 'ARG': 374151, 'TRP': 133431, 'GLY': 511323, 'PHE': 358416, 'ILE': 437133, 'VAL': 495804, 'TYR': 311241, 'SER': 561864, 'LYS': 455313, 'MET': 128106, 'CYS': 197997, 'HIS': 166761, 'PRO': 401658, 'ASP': 397155, 'GLN': 354651, 'ASN': 539043}