import sys
from pathlib import Path
import pandas as pd
import pickle

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()

sys.path.append(str(source_location))
sys.path.append(str(scripts_loc))
from utils import get_sabdab_details


def lacks_epitope_atoms(buried_fullab, pdb_idcode):

    n_atoms_per_chain = [
        len(each) == 0
        for each in buried_fullab[buried_fullab.idcode
                                  == pdb_idcode].epitope_atoms]
    return all(n_atoms_per_chain)


def get_df_dataset(home_dir):
    # df_sabdab_all = pd.read_csv(Path.joinpath(
    #     source_location, 'structures/sabdab_summary_all.tsv'), sep="\t")
    df_sabdab_90 = pd.read_csv(Path.joinpath(
        casa_dir,
        'structures/sabdab_summary_90_May2022.tsv'),
        sep="\t")

    df_buried = pd.read_pickle(Path.joinpath(casa_dir,
                                             'data/epitope_buried.pickle'))

    # df_interactions = pd.read_pickle(Path.joinpath(source_location,
    #                                                'data/interactions.pickle'))

    protein_antigens = df_sabdab_90.query(
        "antigen_type == antigen_type and antigen_type.str.contains('protein')",
        engine='python').drop_duplicates()
    ab_protein_antigens = set(protein_antigens.pdb.values)
    all_saddab_proteins = set(df_sabdab_90.pdb.values)
    print(
        f"SabDab protein antigen:\n"
        f"{len(ab_protein_antigens)} proteins out of {len(all_saddab_proteins)}, "
        f"{round(len(ab_protein_antigens) / len(all_saddab_proteins) * 100, 1)}%")

    ab_both_chains = set(protein_antigens.query(
        "Hchain == Hchain and Lchain == Lchain").pdb.values)
    ab_single_H_chain = set(protein_antigens.query(
        "Hchain == Hchain").pdb.values)
    ab_single_L_chain = set(protein_antigens.query(
        "Lchain == Lchain").pdb.values)

    n_ab_no_Hchain = len(ab_protein_antigens) - len(ab_single_H_chain)
    n_ab_no_Lchain = len(ab_protein_antigens) - len(ab_single_L_chain)

    print(f"All: {len(ab_protein_antigens)}\nNo Hchain: {n_ab_no_Hchain}\nNo Lchain: {n_ab_no_Lchain}\nBoth chains: {len(ab_both_chains)}")

    buried_fullab = df_buried[df_buried.idcode.isin(ab_both_chains)]
    print(
        f"Buried surfaces of {len(set(df_buried.idcode.values))} proteins\n"
        f"with both chains: {len(set(buried_fullab.idcode.values))}"
    )

    """
    # Repeating the epitope_residues column, but as a tuple instead of a list
    buried_fullab.loc[:, 'epitope_residues_tuple'] = buried_fullab.apply(
        lambda row: tuple(row.epitope_residues), axis=1)

    # Now epitope_residues_tuple is hasheable and can be used to discard duplicates
    df_dataset = buried_fullab.drop_duplicates(
        subset=['idcode', 'chain_type', 'cdr', 'cdr_seq',
                'epitope_residues_tuple'])
    """

    return buried_fullab
