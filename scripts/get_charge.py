import numpy as np
import pandas as pd
import os
from pypka import Titration
import traceback

TITRABLE_RESIDUES = ("GLU", "ASP", "HIS", "TYR", "CYS", "LYS")


def get_interface_charge(tit_res, chains, args, tits):
    charge = {chain: np.zeros(len(tits["pHs"])) for chain in chains}
    for res in tit_res:
        chain, resnumb, resname = res

        if resname in TITRABLE_RESIDUES:
            res = (chain, str(resnumb), resname)
            if res in tits:
                charge[chain] += tits[res]
            elif resname != "CYS":
                print("MISSING", chain, resnumb, resname)
                exit()

    for arg in args:
        arg_chain, _ = arg
        charge[arg_chain] += 1

    ph_i = tits["pHs"].index(7.2)
    physiological_charge = 0
    total_charge = np.zeros(len(tits["pHs"]))
    for chain in chains:
        physiological_charge += charge[chain][ph_i]
        total_charge += charge[chain][ph_i]
    physiological_charge = round(physiological_charge, 1)

    return total_charge, physiological_charge


def get_titrable(col):
    titrable = set()
    all_chains = set()
    args = set()
    for residues in df_buried[df_buried.idcode == idcode][col].values:
        for res in residues:
            chain, resnumb, resname, _, _ = res
            if resname in TITRABLE_RESIDUES:
                titrable.add((chain, resnumb, resname))
            elif resname == "ARG":
                args.add((chain, resnumb))
            all_chains.add(chain)

    return titrable, all_chains, args


def titrate(idcode, interface_sites, interface_chains, label):
    pdbfname = f"../structures/raw/{idcode}.pdb"
    newpdbfname = f"../structures/pkas/{idcode}.pdb"
    pkasfname = f"../structures/pkas/{idcode}_tits_{label}.txt"

    if os.path.isfile(pkasfname):
        return read_tit(pkasfname)

    params = {
        "structure": newpdbfname,
        "ncpus": 8,
        "epsin": 15,
        "ionicstr": 0.1,
        "pbc_dimensions": 0,
        "pH": "6.0,8.0",
        "pHstep": 0.2,
        "titration_output": pkasfname,
    }

    sites = {}
    for site in interface_sites:
        chain, resnumb, _ = site
        if chain not in sites:
            sites[chain] = []
        sites[chain].append(resnumb)

    new_pdb_content = ""
    with open(pdbfname) as f:
        for line in f:
            if not line.startswith("ATOM "):
                continue
            atom_chain = line[21]
            if atom_chain in interface_chains:
                new_pdb_content += line
    with open(newpdbfname, "w") as fnew:
        fnew.write(new_pdb_content)

    try:
        Titration(params, sites=sites)
    except Exception:
        with open(pkasfname, "w") as ferr:
            traceback.print_exc(file=ferr)
        return False

    return read_tit(pkasfname)


def read_tit(fname):
    tits = {}
    tmp_tits = []
    with open(fname) as f:
        content = f.readlines()

        for i in content[0].split():
            tmp_tits.append([])

        for line in content[1:]:
            for i, col in enumerate(line.split()):
                tmp_tits[i].append(float(col))

        tits["pHs"] = tmp_tits[0]
        for i, res in enumerate(content[0].split()[1:]):
            resnumb, chain, resname = res.split("_")

            res_charge = np.array(tmp_tits[i + 1])
            if resname in ("ASP", "GLU", "TYR", "CYS"):
                res_charge -= 1
            tits[(chain, resnumb, resname)] = res_charge

    return tits


if __name__ == "__main__":
    df_buried = pd.read_pickle("../data/epitope_buried.pickle")

    idcodes = set(df_buried.idcode.values)

    epi_charges_bound = []
    para_charges_bound = []
    epi_charges_unbound = []
    para_charges_unbound = []
    for idcode in idcodes:
        titrable_epitope, epitope_chains, args_epitope = get_titrable(
            "ag_ab_interface_res"
        )
        titrable_paratope, paratope_chains, args_paratope = get_titrable(
            "ab_ag_interface_res"
        )

        print(idcode, len(titrable_epitope), len(titrable_paratope))

        epitope_tits = titrate(idcode, titrable_epitope, epitope_chains, "epi")
        paratope_tits = titrate(idcode, titrable_paratope, paratope_chains, "para")

        if not epitope_tits or not paratope_tits:
            print("PypKa run failed")
            continue

        epitope_charge, epitope_physiological_charge = get_interface_charge(
            titrable_epitope, epitope_chains, args_epitope, epitope_tits
        )

        paratope_charge, paratope_physiological_charge = get_interface_charge(
            titrable_paratope, paratope_chains, args_paratope, paratope_tits
        )

        print("epitope (unbound)", paratope_physiological_charge)
        print("paratope (unbound)", paratope_physiological_charge)

        epi_charges_unbound.append(epitope_charge)
        para_charges_unbound.append(paratope_charge)

        interface_sites = set.union(titrable_epitope, titrable_paratope)
        interface_chains = set.union(epitope_chains, paratope_chains)
        tits = titrate(idcode, interface_sites, interface_chains, "complex")

        epitope_charge, epitope_physiological_charge = get_interface_charge(
            titrable_epitope, epitope_chains, args_epitope, tits
        )

        paratope_charge, paratope_physiological_charge = get_interface_charge(
            titrable_paratope, paratope_chains, args_paratope, tits
        )

        print("epitope (bound)", paratope_physiological_charge)
        print("paratope (bound)", paratope_physiological_charge)

        epi_charges_bound.append(epitope_charge)
        para_charges_bound.append(paratope_charge)

    

        
    
