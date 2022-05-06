import pandas as pd
import os
import traceback
import pickle
from pypka import Titration
import numpy as np
from scipy.spatial import distance
import plotly.graph_objects as go

from utils import read_pdb_line

TITRABLE_RESIDUES = ("GLU", "ASP", "HIS", "TYR", "CYS", "LYS")
PKAS_MOD = {
    "ASP": 3.94,
    "GLU": 4.25,
    "HIS": 6.54,
    "CYS": 8.55,
    "LYS": 10.40,
    "TYR": 9.84,
    "NTR": 8.00,
    "CTR": 3.67,
}
HBONDS_SITE = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "HIS": ["ND1", "NE2"],
    "CYS": ["SG"],
    "LYS": ["NZ"],
    "TYR": ["OH"],
    "ARG": ["NH2"],
}


def get_titrable(col):
    titrable = set()
    all_chains = set()
    for residues in df_buried[df_buried.idcode == idcode][col].values:
        for res in residues:
            chain, resnumb, resname, _, _ = res
            if resname in TITRABLE_RESIDUES:
                titrable.add((chain, resnumb, resname))
            all_chains.add(chain)
    return titrable, all_chains


def titrate(idcode, interface_sites, interface_chains):
    pdbfname = f"../structures/raw/{idcode}.pdb"
    newpdbfname = f"../structures/pkas/{idcode}.pdb"
    pkasfname = f"../structures/pkas/{idcode}_pkas.txt"

    if os.path.isfile(newpdbfname):
        if os.path.isfile(pkasfname):
            return read_pkas(pkasfname)
        return {}

    params = {
        "structure": newpdbfname,
        "ncpus": 32,
        "epsin": 15,
        "ionicstr": 0.1,
        "pbc_dimensions": 0,
        "output": pkasfname,
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

    return read_pkas(pkasfname)


def read_pkas(fname):
    pkas = {}
    with open(fname) as f:
        for line in f:
            if len(line.split()) != 4:
                break
            resnumb, resname, pka, chain = line.split()
            try:
                pka = float(pka)
            except:
                pka = None
            pkas[(chain, resnumb, resname)] = pka
    return pkas


def calc_diff(pkas):
    diffs = {}
    for site, pka in pkas.items():
        if pka:
            diffs[site] = pka - PKAS_MOD[site[2]]

    return diffs


def average(i):
    return sum(i) / len(i)


def diffs_hists(diffs_epitope, diffs_paratope):
    def sort_by_restype(diffs):
        res_diffs = {}
        for res, diff in diffs.items():
            _, _, _, resname = res
            if resname not in res_diffs:
                res_diffs[resname] = []
            res_diffs[resname].append(diff)
        return res_diffs

    print(len(diffs_epitope), len(diffs_paratope))

    res_diffs_epitope = sort_by_restype(diffs_epitope)
    res_diffs_paratope = sort_by_restype(diffs_paratope)

    for resname in res_diffs_epitope.keys():
        if resname not in res_diffs_paratope or resname not in res_diffs_epitope:
            continue
        print(
            resname,
            average(res_diffs_epitope[resname]),
            average(res_diffs_paratope[resname]),
            len(res_diffs_epitope[resname]),
            len(res_diffs_paratope[resname]),
        )
        res_hist_epitope = np.histogram(
            res_diffs_epitope[resname], range=(-7, 7), bins=13, density=True
        )
        res_hist_paratope = np.histogram(
            res_diffs_paratope[resname], range=(-7, 7), bins=13, density=True
        )
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                y=res_hist_epitope[0],
                x=res_hist_epitope[1] + 0.5,
                mode="lines+markers",
                line_shape="spline",
                name="Epitope",
            )
        )
        fig.add_trace(
            go.Scatter(
                y=res_hist_paratope[0],
                x=res_hist_paratope[1] + 0.5,
                mode="lines+markers",
                line_shape="spline",
                name="Paratope",
            )
        )
        fig.update_layout(
            xaxis_title="p<i>K</i><sub>a</sub> Shift",
            yaxis_title="Probability Density",
            template="plotly_white",
            xaxis_range=[-7.1, 7.1],
            autosize=False,
            width=500,
            height=500,
            margin=dict(l=50, r=50, b=100, t=100, pad=4),
        )
        fig.write_image(f"../plots/pkas/{resname}_hist.png")


def get_structure_distance(idcode, titrable_epitope, titrable_paratope):
    atoms = []
    with open(f"../structures/pkas/{idcode}.pdb") as f:
        for line in f:
            if line.startswith("ATOM "):
                (
                    chain,
                    cdr_id,
                    (resname, resnumb),
                    (aname, anumb),
                    coord,
                ) = read_pdb_line(line)
                if cdr_id:
                    resnumb = f"{resnumb}{cdr_id}"
                if aname[0] in "NOS":
                    atom = (chain, resname, resnumb, aname, anumb, coord)
                    atoms.append(atom)

    df_atoms = pd.DataFrame(
        atoms, columns=("chain", "resname", "resnumb", "aname", "anumb", "coord")
    )

    dist_matrix = distance.cdist(
        df_atoms["coord"].values.tolist(),
        df_atoms["coord"].values.tolist(),
    )

    wm_range = dist_matrix < 3.5
    hbonds = {}

    for i1, a1 in enumerate(wm_range):
        if a1.sum() == 0:
            continue
        for i2 in a1.nonzero()[0]:
            atom1 = df_atoms.iloc[i1]
            atom2 = df_atoms.iloc[i2]
            if atom1.resnumb != atom2.resnumb:
                chain1, resnumb1, resname1 = (
                    atom1.chain,
                    atom1.resnumb,
                    atom1.resname,
                )
                chain2, resnumb2, resname2 = (
                    atom2.chain,
                    atom2.resnumb,
                    atom2.resname,
                )

                katom1 = (chain1, resnumb1, resname1)
                katom2 = (chain2, resnumb2, resname2)

                if (
                    katom1 in titrable_epitope
                    or katom1 in titrable_paratope
                    and atom1.aname in HBONDS_SITE[resname1]
                ):
                    interaction = (
                        katom2,
                        atom1.aname,
                        atom2.aname,
                        round(dist_matrix[i1][i2], 1),
                    )
                    if katom1 not in hbonds:
                        hbonds[katom1] = []
                    hbonds[katom1].append(interaction)
                if (
                    katom2 in titrable_epitope
                    or katom2 in titrable_paratope
                    and atom2.aname in HBONDS_SITE[resname2]
                ):
                    interaction = (
                        katom1,
                        atom2.aname,
                        atom1.aname,
                        round(dist_matrix[i1][i2], 1),
                    )
                    if katom2 not in hbonds:
                        hbonds[katom2] = []
                    hbonds[katom2].append(interaction)

    return hbonds


def hbonds_anal(hbonds, titrable_epitope, titrable_paratope):
    # from pprint import pprint

    salt_bridges = [0, 0, 0, 0, 0, 0]
    para_salt_bridges_byres = {
        "GLU": 0,
        "ASP": 0,
        "HIS": 0,
        "TYR": 0,
        "CYS": 0,
        "LYS": 0,
    }
    epi_salt_bridges_byres = {
        "GLU": 0,
        "ASP": 0,
        "HIS": 0,
        "TYR": 0,
        "CYS": 0,
        "LYS": 0,
    }

    from pprint import pprint

    # WRONG WAY OF COUNTING SALT BRIDGES!
    # ONLY res2 type is being considered

    for res in titrable_paratope:
        # print(res)
        if res in hbonds:
            res_hbonds = hbonds[res]
            # pprint(res_hbonds)
            trigger_saltb = False
            for hbond in res_hbonds:
                (chain, _, resname2), aname1, aname2, _ = hbond

                if (
                    chain in epitope_chains
                    and resname2 in HBONDS_SITE
                    and aname2 in HBONDS_SITE[resname2]
                ):
                    salt_bridges[2] += 1
                    para_salt_bridges_byres[res[2]] += 1
                    trigger_saltb = True
                    break
            if not trigger_saltb:
                salt_bridges[1] += 1
        else:
            salt_bridges[0] += 1
        # print(salt_bridges)

    for res in titrable_epitope:
        # print(res)
        if res in hbonds:
            res_hbonds = hbonds[res]
            # pprint(res_hbonds)
            trigger_saltb = False
            for hbond in res_hbonds:
                (chain, _, resname2), aname1, aname2, _ = hbond
                if (
                    chain in paratope_chains
                    and resname2 in HBONDS_SITE
                    and aname2 in HBONDS_SITE[resname2]
                ):
                    salt_bridges[5] += 1
                    epi_salt_bridges_byres[res[2]] += 1
                    trigger_saltb = True
                    break
            if trigger_saltb:
                salt_bridges[4] += 1
        else:
            salt_bridges[3] += 1

    print(salt_bridges)
    print(para_salt_bridges_byres)
    print(epi_salt_bridges_byres)
    return salt_bridges, para_salt_bridges_byres, epi_salt_bridges_byres


if __name__ == "__main__":
    df_buried = pd.read_pickle("../data/epitope_buried.pickle")

    idcodes = []
    with open("../final_ids.txt") as f:
        for line in f:
            idcodes.append(line.strip())
    print(len(idcodes))

    diffs_epitope = {}
    diffs_paratope = {}
    salt_bridges, para_salt_bridges_byres, epi_salt_bridges_byres = [], [], []
    for idcode in idcodes:
        titrable_epitope, epitope_chains = get_titrable("ag_ab_interface_res")
        titrable_paratope, paratope_chains = get_titrable("ab_ag_interface_res")

        print(idcode, len(titrable_epitope), len(titrable_paratope))

        interface_sites = set.union(titrable_epitope, titrable_paratope)
        interface_chains = set.union(epitope_chains, paratope_chains)

        pkas = titrate(idcode, interface_sites, interface_chains)

        prot_diffs = calc_diff(pkas)
        print([round(i, 1) for i in prot_diffs.values()])

        for res, diff in prot_diffs.items():
            resid = (idcode, res[0], res[1], res[2])
            if res[0] in epitope_chains:
                diffs_epitope[resid] = diff
            elif res[0] in paratope_chains:
                diffs_paratope[resid] = diff
            else:
                raise Exception

        if pkas:
            prot_hbonds = get_structure_distance(
                idcode, titrable_epitope, titrable_paratope
            )

            (
                prot_salt_bridges,
                prot_para_salt_bridges_byres,
                prot_epi_salt_bridges_byres,
            ) = hbonds_anal(prot_hbonds, titrable_epitope, titrable_paratope)
            salt_bridges.append(prot_salt_bridges)
            para_salt_bridges_byres.append(prot_para_salt_bridges_byres)
            epi_salt_bridges_byres.append(prot_epi_salt_bridges_byres)

    # diffs_hists(diffs_epitope, diffs_paratope)
