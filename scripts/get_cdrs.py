import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1 as AA_CONVERTER
from abnumber import Chain
from utils import read_pdb_line, get_sabdab_details, remove_redundant
import numpy as np
from scipy.spatial import distance_matrix
from multiprocessing import Pool


class Structure:
    def __init__(self, idcode, df_pdb):
        self.idcode = idcode
        self.pdb_path = f"../structures/raw/{idcode}.pdb"

        self.chains = []
        for i, values in df_pdb.iterrows():
            chains = list(values[["Hchain", "Lchain"]].dropna().unique())
            if len(chains) == 2:
                self.chains += chains

        self.antigen_chains = df_pdb.antigen_chain.values[0].replace(" | ", "")

        self.get_sequence()
        self.ab_chains = {}

    def get_sequence(self):
        self.sequence = {chain: "" for chain in self.chains}
        self.seqnumbs = {chain: [] for chain in self.chains}
        last = None
        with open(self.pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain, res, *_ = read_pdb_line(line)
                    (resname, resnumb) = res
                    new = resname, chain, resnumb
                    if chain in self.chains and last != new:
                        res = AA_CONVERTER[resname.capitalize()]
                        self.sequence[chain] += res
                        self.seqnumbs[chain].append(resnumb)
                    last = new[:]

    def get_cdrs(self, scheme):
        cdrs = {}
        for ichain in range(0, len(self.chains), 2):
            hchain, lchain = self.chains[ichain : ichain + 2]
            chains = {}
            for chainid in (hchain, lchain):
                seq = self.sequence[chainid]
                print(chainid, seq, scheme)
                try:
                    chain = Chain(seq, scheme=scheme)
                    chains[chainid] = chain
                except:
                    continue

            if len(chains) != 2:
                continue

            for chainid in (hchain, lchain):
                seq = self.sequence[chainid]
                chain = chains[chainid]
                ab_chain = AbChain(seq, chainid, chain.chain_type)
                self.ab_chains[chainid] = ab_chain

                cdrs[chainid] = {
                    1: chain.cdr1_seq,
                    2: chain.cdr2_seq,
                    3: chain.cdr3_seq,
                }
                for i, cdr in cdrs[chainid].items():
                    b = seq.index(cdr)
                    e = b + len(cdr) - 1

                    b_resnumb = self.seqnumbs[chainid][b]
                    e_resnumb = self.seqnumbs[chainid][e]

                    ab_chain.add_cdr(i, cdr, b_resnumb, e_resnumb)

    def get_cdr_epitope(self, cutoff):
        if not self.ab_chains:
            return

        cur_cdr = None
        cdr_i = None
        prev_resnumb = None
        # Assign atoms to CDRs
        with open(self.pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain, res, atom, coords = read_pdb_line(line)
                    if chain in self.ab_chains:
                        (resname, resnumb) = res

                        abchain = self.ab_chains[chain]

                        if not cur_cdr:
                            for i, cdr in abchain.cdrs.items():
                                if resnumb >= cdr.i_begin and resnumb <= cdr.i_end:
                                    if resnumb == cdr.i_begin:
                                        cur_cdr = i
                                        cdr_i = -1
                                        break

                        if cur_cdr:
                            cdr = abchain.cdrs[cur_cdr]

                            if prev_resnumb != resnumb:
                                cdr_i += 1

                            if resnumb > cdr.i_end or cdr_i >= len(cdr.seq):
                                # print("Reset")
                                prev_resnumb = None
                                cur_cdr = None
                            else:
                                # print(
                                #     resnumb,
                                #     cdr.i_begin,
                                #     cdr.i_end,
                                #     resname,
                                #     AA_CONVERTER[resname.capitalize()],
                                #     cdr.seq[cdr_i],
                                #     cdr_i,
                                # )
                                assert (
                                    AA_CONVERTER[resname.capitalize()] == cdr.seq[cdr_i]
                                ), f"CDR from sequence does not match the structure\n{AA_CONVERTER[resname.capitalize()]}\n{cdr.seq[cdr_i]}"
                                cdr.add_atom(atom, res, coords)

                                prev_resnumb = resnumb

        candidates_details = []
        candidates_coords = []
        with open(self.pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain, res, atom, coords = read_pdb_line(line)
                    if chain in self.antigen_chains:
                        (resname, resnumb) = res

                        candidates_details.append((chain, res, atom))
                        candidates_coords.append(coords)

        if len(candidates_coords) == 0:
            return

        residues = []
        for chainid, abchain in self.ab_chains.items():
            for cdr in abchain.cdrs.values():
                cdr_atoms = [atom.coords for atom in cdr.atoms.values()]

                to_print = {}
                for chain, res, atom in candidates_details:
                    if chain not in to_print:
                        to_print[chain] = []
                    to_print[chain].append(atom[1])

                if len(cdr_atoms):
                    dm = distance_matrix(candidates_coords, cdr_atoms)
                else:
                    continue
                cutoff_dm = dm < cutoff
                for i in range(dm.shape[0]):
                    if cutoff_dm[i].sum() > 0:
                        anumb = candidates_details[i][2][1]
                        resnumb = candidates_details[i][1][1]
                        cdr.save_epitope(anumb, resnumb)

                # cdr.print_pymol_selection()
                # exit()

    def print_cdrs(self):
        for chain in self.ab_chains.values():
            for i, cdr in chain.iter_cdrs():
                print(
                    chain.chain_id, chain.chain_type, i, cdr.seq, cdr.i_begin, cdr.i_end
                )

    def save_cdr_epitopes(self):
        cdrs = []
        for chain in self.ab_chains.values():
            for i, cdr in chain.iter_cdrs():
                cdrs.append(
                    [
                        self.idcode,
                        chain.chain_id,
                        chain.chain_type,
                        i,
                        cdr.seq,
                        cdr.i_begin,
                        cdr.i_end,
                        [atom.anumb for atom in cdr.atoms.values()],
                        cdr.epitope,
                        cdr.epitope_res,
                    ]
                )

        return cdrs


class AbChain:
    def __init__(self, seq, chain_id, chain_type):
        self.seq = seq
        self.chain_id = chain_id
        self.chain_type = chain_type

        self.cdrs = {}

    def add_cdr(self, cdr_type, seq, begin, end):
        cdr = CDR(cdr_type, seq, begin, end)
        self.cdrs[cdr_type] = cdr

    def iter_cdrs(self):
        for i, cdr in self.cdrs.items():
            yield i, cdr


class CDR:
    def __init__(self, cdr_type, seq, b, e):
        self.numb = cdr_type
        self.seq = seq
        self.i_begin = b
        self.i_end = e
        self.size = len(seq)
        self.atoms = {}
        self.epitope = []
        self.epitope_res = []

    def add_atom(self, atom, res, coords):
        new_atom = Atom(atom, res, coords)
        self.atoms[atom[1]] = new_atom

    def save_epitope(self, anumb, resnumb):
        self.epitope.append(anumb)
        self.epitope_res.append(resnumb)

    def print_pymol_selection(self):
        sel_cdr = "+".join([str(i) for i in self.atoms.keys()])
        sel_epitope = "+".join([str(i) for i in self.epitope])
        print("CDR:", len(self.atoms), sel_cdr)
        print("Epitope:", len(self.epitope), sel_epitope)


class Atom:
    def __init__(self, atom, res, coords):
        self.resname, self.resnumb = res
        self.aname, self.anumb = atom
        self.coords = np.array(coords)


def run(x):
    i, pdb = x
    print(i, len(unique_pdbs), pdb)
    df_pdb = df_details.query(f'pdb == "{pdb}"')[["Hchain", "Lchain", "antigen_chain"]]

    # try:
    struct = Structure(pdb, df_pdb)
    struct.get_cdrs(scheme)
    struct.print_cdrs()
    struct.get_cdr_epitope(cutoff)
    cdrs = struct.save_cdr_epitopes()

    # except:
    #     cdrs = None
    #     print("failed ->", pdb)

    if cdrs:
        df = pd.DataFrame(
            cdrs,
            columns=(
                "idcode",
                "chainID",
                "chain_type",
                "cdr",
                "cdr_seq",
                "cdr_begin",
                "cdr_end",
                "cdr_atoms",
                "epitope_atoms",
                "epitope_residues",
            ),
        )
        df, removed = remove_redundant(df)
        return df, removed
    else:
        print("skipped")


if __name__ == "__main__":
    df_details = get_sabdab_details()
    scheme = "chothia"
    cutoff = 5

    df = pd.DataFrame(
        columns=(
            "idcode",
            "chainID",
            "chain_type",
            "cdr",
            "cdr_seq",
            "cdr_begin",
            "cdr_end",
            "cdr_atoms",
            "epitope_atoms",
            "epitope_residues",
        ),
    )
    df_removed = df.copy()
    unique_pdbs = df_details.pdb.unique()

    # df_r = run((0, "1adq"))
    # print(
    #     df_r[0][
    #         [
    #             "idcode",
    #             "chainID",
    #             "chain_type",
    #             "cdr",
    #             "cdr_atoms",
    #             "epitope_residues",
    #         ]
    #     ]
    # )
    # exit()

    # for i, pdb in enumerate(unique_pdbs[675:]):
    #     print(f"### {i} ###")
    #     df_r = run((i, pdb))
    #     print(df_r[["idcode", "chainID", "chain_type", "cdr", "epitope_residues"]])
    #     exit()

    with Pool(8) as p:
        rs = p.map(run, [(i, pdb) for i, pdb in enumerate(unique_pdbs)])

    for r in rs:
        if r[0] is not None:
            df = df.append(r[0], ignore_index=True)
        if r[1] is not None:
            df_removed = df_removed.append(r[1], ignore_index=True)

    print(f"PDBs: {len(df.idcode.unique())}\tCDRs: {len(df)}\tAbAg: {len(df)/6}")
    df.to_pickle("../data/cdr_epitope.pickle")

    print(
        f"PDBs: {len(df_removed.idcode.unique())}\tCDRs: {len(df_removed)}\tAbAg: {len(df_removed)/6}"
    )
    df_removed.to_pickle("../data/cdr_epitope_redundant.pickle")
