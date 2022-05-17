import pickle
import sys
from pathlib import Path
import glob
import numpy as np
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])

source_location = (Path(__file__) / '..').resolve()
scripts_loc = (Path(__file__) / '..' / '..').resolve()
sys.path.append(source_location)

casa_dir = (Path(__file__) / '..' / '..' / '..').resolve()
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")


def get_chain_from_biggest_interface(df_pdb):
    top_interface = 0
    if len(df_pdb.chainID) > 6:
        interface_sizes_cdr = []
        for row in df_pdb.iterrows():
            interface_sizes_cdr.append(
                len(row[1].ag_ab_interface_res) + len(row[1].ab_ag_interface_res))
        assert len(interface_sizes_cdr) % 6 == 0, f"BAD: {pdb_idcode}. "\
            "There should be a record for each CDR."

        top_interface = interface_sizes_cdr.index(max(interface_sizes_cdr))

    return df_pdb.chainID.values[top_interface]


if __name__ == '__main__':
    print("Reading buried_interface_res.pickle")
    with open(Path.joinpath(
            casa_dir, "data", 'buried_interface_res.pickle'), 'rb') as file:
        buried_interface_res = pickle.load(file)
    print("Reading epitope_buried_cleaned.pickle")
    with open(Path.joinpath(
            casa_dir, "data", 'epitope_buried_cleaned.pickle'), 'rb') as file:
        epitope_buried_cleaned = pickle.load(file)
    print("Starting now.")

    pdb_list = list(set(buried_interface_res.idcode))
    pdb_list.sort()
    filename_dict = {}
    chains_dict = {}
    # check_pdb = '2qad'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        df_pdb = epitope_buried_cleaned.query(f"idcode == '{pdb_idcode}'")
        assert df_pdb.chain_type.values[0] == 'H', f"BAD: {pdb_idcode}. "\
            "First record should correspond to a heavy chain."
        target_heavy_chain = get_chain_from_biggest_interface(df_pdb)

        pdb_string = pdb_idcode + "_complex_??_*.pdb"
        chains_list = glob.glob(str(Path.joinpath(
            casa_dir, "structures", "exposed", pdb_idcode, pdb_string)))
        assert len(chains_list) != 0, f"BAD: {pdb_idcode}. No files. "

        for pdb_full_path in chains_list:
            pdb_filename = Path(pdb_full_path).name
            H_chain = pdb_filename[13]
            L_chain = pdb_filename[14]
            AG_chain_a = pdb_filename[16]
            AG_chain_b = pdb_filename[17]
            # Thank god the antigen can have up to 5 chains,
            # or else I'd have to actually code right
            AG_chain_c = pdb_filename[18] if pdb_filename[18].isupper() else '.'
            AG_chain_d = pdb_filename[19] if pdb_filename[19].isupper() else '.'
            AG_chain_e = pdb_filename[20] if pdb_filename[20].isupper() else '.'
            if H_chain == target_heavy_chain:
                chains_dict[pdb_idcode] = Chains(antibody=(H_chain, L_chain),
                                                 antigen=(
                    AG_chain_a, AG_chain_b, AG_chain_c, AG_chain_d, AG_chain_e))

                filename_dict[pdb_idcode] = pdb_filename
                break
        else:
            raise RuntimeError(f"BAD: {pdb_idcode}. "
                               "Available files are not included in the dataset. ")

    with open(Path.joinpath(data_dir, 'filenames.pkl'), 'wb') as file:
        pickle.dump(filename_dict, file)
    with open(Path.joinpath(data_dir, 'chains.pkl'), 'wb') as file:
        pickle.dump(chains_dict, file)
    with(open(Path.joinpath(data_dir, 'pdb.list'), 'w')) as file:
        all(file.write(f"{pdb}\n") for pdb in pdb_list)
