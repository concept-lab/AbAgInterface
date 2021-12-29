from scipy.spatial import distance
import pandas as pd
import numpy as np
from utils import pqr2xyzr, read_pdb_line, get_sabdab_details

# https://pubs.acs.org/doi/10.1021/ja034729u


class AbAgComplex:
    def __init__(
        self, idcode: str, abh_chain: str, abl_chain: str, ag_chain: str
    ) -> None:
        print("\n", idcode)
        self.idcode = idcode
        self.ab_h_chain = abh_chain
        self.ab_l_chain = abl_chain
        self.ag_chain = ag_chain

        self.ab_h_atoms = None
        self.ab_l_atoms = None
        self.ag_atoms = None

        self.exposed_atoms = []

        # Direct any < 5
        self.direct_ab_h_ag = []
        self.direct_ab_l_ag = []

        # Hbonds NOS < 3.5
        self.hbond_ab_h_ag = []
        self.hbond_ab_l_ag = []

        # Water-mediated NOS < 6-8
        self.wm_ab_h_ag = []
        self.wm_ab_l_ag = []

        self.label_chains()

    def label_chains(self):
        interesting_chains = self.ab_h_chain + self.ab_l_chain + self.ag_chain
        atoms = []
        with open(f"../structures/raw/{self.idcode}.pdb") as f:
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
                    if chain not in interesting_chains:
                        continue
                    atom = (chain, resname, resnumb, aname, anumb, coord)
                    atoms.append(atom)

        df_atoms = pd.DataFrame(
            atoms, columns=("chain", "resname", "resnumb", "aname", "anumb", "coord")
        )
        self.ab_h_atoms = df_atoms.query(f"chain == '{self.ab_h_chain}'")
        self.ab_l_atoms = df_atoms.query(f"chain == '{self.ab_l_chain}'")
        self.ag_atoms = df_atoms.query(f"chain.isin(list('{self.ag_chain}'))")

        print(len(self.ab_h_atoms), len(self.ab_l_atoms), len(self.ag_atoms))

    def find(self):
        print(
            f"Ab H: {self.ab_h_chain}\t\tAb L: {self.ab_l_chain}\t\tAg: {self.ag_chain}"
        )

        # Calculate distances
        abh_ag = (
            distance.cdist(
                self.ab_h_atoms["coord"].values.tolist(),
                self.ag_atoms["coord"].values.tolist(),
            )
            if not self.ab_h_atoms.empty
            else np.array([])
        )
        abl_ag = (
            distance.cdist(
                self.ab_l_atoms["coord"].values.tolist(),
                self.ag_atoms["coord"].values.tolist(),
            )
            if not self.ab_l_atoms.empty
            else np.array([])
        )

        # Water-mediated
        self.get_exposed()
        if not self.exposed_atoms:
            return True

        self.wm_ab_h_ag = self.find_watermed(abh_ag, self.ab_h_atoms)
        self.wm_ab_l_ag = self.find_watermed(abl_ag, self.ab_l_atoms)
        print(
            "Water-mediated\t"
            f"{len(self.wm_ab_h_ag['Ab_i'].unique()):5}|{len(self.wm_ab_l_ag['Ab_i'].unique()):5>}\t"
            f"{len(self.wm_ab_h_ag['Ag_i'].unique()):5}|{len(self.wm_ab_l_ag['Ag_i'].unique()):5>}\t"
            f"{len(self.wm_ab_h_ag):5}|{len(self.wm_ab_l_ag)}"
        )

        # H bonds
        self.hbond_ab_h_ag = self.find_hbonds(abh_ag, self.ab_h_atoms)
        self.hbond_ab_l_ag = self.find_hbonds(abl_ag, self.ab_l_atoms)
        print(f"Hbonds\t\t{len(self.hbond_ab_h_ag):5}|{len(self.hbond_ab_l_ag)}")

        # Direct Interactions
        self.direct_ab_h_ag = self.find_direct(abh_ag)
        self.direct_ab_l_ag = self.find_direct(abl_ag)
        print(f"Direct\t\t{len(self.direct_ab_h_ag):5}|{len(self.direct_ab_l_ag)}")

        return False

    def find_direct(self, ab_ag_dists):
        wm_range = ab_ag_dists < 5
        direct = []
        for i_ab, ab_atom in enumerate(wm_range):
            if ab_atom.sum() > 0:
                for i_ag in ab_atom.nonzero()[0]:
                    interaction = (
                        i_ab,
                        i_ag,
                        round(ab_ag_dists[i_ab][i_ag], 1),
                    )
                    direct.append(interaction)

        return pd.DataFrame(direct, columns=("Ab_i", "Ag_i", "distance"))

    def find_hbonds(self, ab_ag_dists, ab_atoms):
        wm_range = ab_ag_dists < 3.5
        hbonds = []
        for i_ab, ab_atom in enumerate(wm_range):
            if ab_atom.sum() > 0 and ab_atoms.iloc[i_ab]["aname"][0] in "NOS":
                for i_ag in ab_atom.nonzero()[0]:
                    if self.ag_atoms.iloc[i_ag]["aname"][0] in "NOS":
                        interaction = (
                            i_ab,
                            i_ag,
                            round(ab_ag_dists[i_ab][i_ag], 1),
                        )
                        hbonds.append(interaction)

        return pd.DataFrame(hbonds, columns=("Ab_i", "Ag_i", "distance"))

    def find_watermed(self, ab_ag_dist, ab_atoms):
        wm_range = (ab_ag_dist > 6) & (ab_ag_dist < 8)
        interactions = []
        for i_atom, ab_atom in enumerate(wm_range):
            if ab_atom.sum() > 0 and ab_atoms.iloc[i_atom]["aname"][0] in "NOS":
                ab_atom_res = ab_atoms.iloc[i_atom][
                    ["chain", "aname", "resnumb"]
                ].values.tolist()

                if ab_atom_res not in self.exposed_atoms:
                    continue

                for i_ag in ab_atom.nonzero()[0]:
                    if self.ag_atoms.iloc[i_ag]["aname"][0] in "NOS":
                        ag_atom_res = self.ag_atoms.iloc[i_ag][
                            ["chain", "aname", "resnumb"]
                        ].values.tolist()
                        if ag_atom_res in self.exposed_atoms:
                            interaction = (
                                i_atom,
                                i_ag,
                                round(ab_ag_dist[i_atom][i_ag], 1),
                            )
                            interactions.append(interaction)

        return pd.DataFrame(interactions, columns=("Ab_i", "Ag_i", "distance"))

    def get_exposed(self):
        from glob import glob

        pqr_path = f"../structures/exposed/{self.idcode}/{self.idcode}_complex_*.pqr"
        pqr = glob(pqr_path)
        if len(pqr) != 1:
            return None

        exposed_path = (
            f"../structures/exposed/{self.idcode}/{self.idcode}_complex_*.txt"
        )
        exposed = glob(exposed_path)
        if len(exposed) != 1:
            return None

        atoms = []
        with open(pqr[0]) as f:
            for line in f:
                if line.startswith("ATOM "):
                    chain, cdr_id, (_, resnumb), (aname, anumb), coord = read_pdb_line(
                        line
                    )
                    if cdr_id:
                        resnumb = f"{resnumb}{cdr_id}"
                    atoms.append([chain, aname, resnumb])

        with open(exposed[0]) as f:
            for line in f:
                i_atom = int(line)
                atom = atoms[i_atom]
                self.exposed_atoms.append(atom)

    def save(self):
        return (
            self.idcode,
            self.ab_h_chain,
            self.ab_l_chain,
            self.ag_chain,
            [
                self.ab_h_atoms,
                self.ab_l_atoms,
                self.ag_atoms,
                self.exposed_atoms,
                self.direct_ab_h_ag,
                self.direct_ab_l_ag,
                self.hbond_ab_h_ag,
                self.hbond_ab_l_ag,
                self.wm_ab_h_ag,
                self.wm_ab_l_ag,
            ],
        )


def get_interactions(fname):
    df_sabdab_all = get_sabdab_details()

    df = pd.read_pickle(fname)
    interactions = []
    idcodes = sorted(set(df["idcode"].values))
    for i_idcode, idcode in enumerate(idcodes):
        print(i_idcode, idcode)
        sabdab_prot = df_sabdab_all.query(f'pdb == "{idcode}"')
        antigen_type = sabdab_prot["antigen_type"].values[0]
        if "peptide" in antigen_type:
            continue
        hchain = (
            "".join([str(i) if i == i else "" for i in sabdab_prot["Hchain"].values])
            if sabdab_prot["Hchain"].values.any()
            else ""
        )
        lchain = (
            "".join([str(i) if i == i else "" for i in sabdab_prot["Lchain"].values])
            if sabdab_prot["Lchain"].values.any()
            else ""
        )
        ag_chain = "".join(
            [i.replace(" | ", "") for i in sabdab_prot["antigen_chain"].values]
        )

        wm = AbAgComplex(idcode, hchain, lchain, ag_chain)
        skip = wm.find()
        if skip:
            continue
        data = wm.save()
        interactions.append(data)

    import pickle

    with open("../data/interactions.pickle", "wb") as handle:
        pickle.dump(interactions, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    f_epitopes = "../data/cdr_epitope.pickle"
    get_interactions(f_epitopes)
