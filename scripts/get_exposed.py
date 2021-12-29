import os
import pandas as pd
from utils import read_pdb_line, get_sabdab_details, pqr2xyzr
import subprocess

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
EXPOSED_DIR = f"{CUR_DIR}/../structures/exposed/"


def iterate_prots(df_fname):
    df = pd.read_pickle(df_fname)  # .iloc[::-1]
    for i, row in df.iterrows():
        yield row
    # exit()


class System:
    def __init__(self, fname, atoms):
        self.fname = fname
        self.atoms = atoms
        self.indices = {}
        self.surface_atoms = []
        self.fpdb = f"{PROT_DIR}{fname}.pdb"
        self.fpqr = f"{PROT_DIR}{fname}.pqr"
        self.fxyzr = f"{PROT_DIR}{fname}.xyzr"
        self.fexposed = f"{PROT_DIR}{fname}_exposed.txt"

    def trim_xyzr(self, complex_pqr, complex_xyzr):
        self_xyzr = []
        idx = 0
        with open(complex_pqr) as f_pqr, open(complex_xyzr) as f_xyzr:
            for line_pqr, line_xyzr in zip(f_pqr, f_xyzr):
                chain, _, (_, resnumb), (_, anumb), (x, y, z), radius = read_pdb_line(
                    line_pqr, pqr=True
                )
                if chain in self.atoms and resnumb in self.atoms[chain]:
                    self_xyzr.append(line_xyzr)
                    self.indices[idx] = anumb
                    idx += 1

        with open(self.fxyzr, "w") as f:
            f.write("".join(self_xyzr))

    def get_surface_atoms(self):
        with open(self.fexposed) as f:
            for line in f:
                idx = line.strip()
                if idx:
                    idx = int(idx)
                    anumb = self.indices[idx]
                    self.surface_atoms.append(anumb)
        self.surface_atoms = set(self.surface_atoms)


def run_command(command):
    try:
        cmd = command
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"-> {command} <- did not run successfully\nMessage: {e.stderr.decode('ascii')}"
        )


def get_systems(idcode, cdr_atoms, antigen_chains, cdr_chain, ab_chains, cdr_numb):
    complex_pdb = []
    ab_cdr_atoms = []
    chain_res = {}
    print(f"{base}/{idcode}.pdb")
    with open(f"{base}/{idcode}.pdb") as f:
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

            if anumb in cdr_atoms and resnumb not in ab_cdr_atoms:
                ab_cdr_atoms.append(resnumb)

            if chain in ab_chains or chain in antigen_chains:
                complex_pdb.append(line)

    if ab_chains[0] != ab_chains[0] or ab_chains[1] != ab_chains[1]:
        otherab_chain = ""
        chains = ab_chains[0] if ab_chains[1] != ab_chains[1] else ab_chains[1]
    else:
        chains = "".join(ab_chains)
        otherab_chain = chains.replace(cdr_chain, "").strip()

    f_complex = f"{idcode}_complex_{chains}_{antigen_chains}"
    f_chaincomplex = f"{idcode}_chaincomplex_{cdr_chain}_{antigen_chains}"
    f_cdrcomplex = f"{idcode}_cdrcomplex_{cdr_chain}{cdr_numb}_{antigen_chains}"
    f_antibody = f"{idcode}_antibody_{chains}"
    f_antigen = f"{idcode}_antigen_{antigen_chains}"

    ab_ag_complex = System(
        f_complex,
        {
            chain: list(set(chain_res[chain]))
            for chain in antigen_chains + cdr_chain + otherab_chain
        },
    )
    cdrchain_ag_complex = System(
        f_chaincomplex,
        {chain: list(set(chain_res[chain])) for chain in antigen_chains + cdr_chain},
    )
    cdr_ag_complex = System(
        f_cdrcomplex,
        dict(
            {chain: list(set(chain_res[chain])) for chain in antigen_chains},
            **{cdr_chain: ab_cdr_atoms},
        ),
    )
    # cdrchain = System(f_cdrchain, {cdr_chain: ab_cdrchain_res})
    antibody = System(
        f_antibody,
        {chain: list(set(chain_res[chain])) for chain in cdr_chain + otherab_chain},
    )
    antigen = System(
        f_antigen, {chain: list(set(chain_res[chain])) for chain in antigen_chains}
    )

    if not os.path.isfile(ab_ag_complex.fpqr):
        with open(ab_ag_complex.fpdb, "w") as f:
            f.write("".join(complex_pdb))

        cmd = f'pdb2pqr30 --ff=AMBER --keep-chain --noopt --nodebump {ab_ag_complex.fpdb} {ab_ag_complex.fpqr}; grep "ATOM" {ab_ag_complex.fpqr} > tmp.pqr; mv tmp.pqr {ab_ag_complex.fpqr}'
        run_command(cmd)

        # check if pdb2pqr crashed
        empty_trigger = False
        with open(ab_ag_complex.fpqr) as f:
            if len(f.readlines()) == 0:
                empty_trigger = True
        if empty_trigger:
            os.system(f"rm -f {ab_ag_complex.fpqr}")
            return [None] * 5
    else:
        print(f"SKIPPING pdb2pqr run of {ab_ag_complex.fpdb}")

    return (ab_ag_complex, cdrchain_ag_complex, cdr_ag_complex, antibody, antigen)


def run_nanoshaper(systems: list):
    for system in systems:
        print(system.fname)
        system.trim_xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr)

        if os.path.isfile(system.fexposed):
            print("SKIPPED", system.fexposed)
            system.get_surface_atoms()
            continue

        cmd = f"""sed "s/__XYZRFILENAME__/{system.fxyzr.split('/')[-1]}/" template.prm > {PROT_DIR}conf.prm"""
        os.system(cmd)

        cmd = f"docker run -i --mount type=bind,source={PROT_DIR},target=/App registry-gitlab.iit.it/sdecherchi/nanoshaper:0.7.8 conf.prm"
        run_command(cmd)

        os.system(
            f"cp {PROT_DIR}exposedIndices.txt {system.fexposed}; rm -f {PROT_DIR}exposedIndices.txt"
        )
        system.get_surface_atoms()
    return


def get_interfaces(atom_res):
    ag_ab_interface = antigen.surface_atoms - ab_ag_complex.surface_atoms
    ag_cdrchain_interface = antigen.surface_atoms - cdrchain_ag_complex.surface_atoms
    ag_cdr_interface = antigen.surface_atoms - cdr_ag_complex.surface_atoms
    ab_ag_interface = antibody.surface_atoms - ab_ag_complex.surface_atoms

    ag_ab_interface_res = [atom_res[i] for i in ag_ab_interface]
    ag_cdrchain_interface_res = [atom_res[i] for i in ag_cdrchain_interface]
    ag_cdr_interface_res = [atom_res[i] for i in ag_cdr_interface]
    ab_ag_interface_res = [atom_res[i] for i in ab_ag_interface]

    # print(set(ag_ab_interface_res))
    # print(set(ag_cdrchain_interface_res))
    # print(set(ag_cdr_interface_res))
    # print(set(ab_ag_interface_res))
    # exit()

    return (
        (ag_ab_interface, ag_ab_interface_res),
        (ag_cdrchain_interface, ag_cdrchain_interface_res),
        (ag_cdr_interface, ag_cdr_interface_res),
        (ab_ag_interface, ab_ag_interface_res),
    )





if __name__ == "__main__":
    df_sabdab_all = get_sabdab_details()

    base = "../structures/raw/"
    f_epitopes = "../data/cdr_epitope.pickle"

    new_df = []
    c = 0
    for prot in iterate_prots(f_epitopes):
        if prot["idcode"] == "6kn9":
            continue

        c += 1
        print(c)

        ab_chain_type = prot["chain_type"]
        if prot["chain_type"] == "K":
            ab_chain_type = "L"

        df_sabdab_prot = df_sabdab_all.query(
            f"pdb == '{prot['idcode']}' and {ab_chain_type}chain == '{prot['chainID']}'"
        )
        antigen_type = df_sabdab_prot["antigen_type"].values[0]
        antigen_chains = df_sabdab_prot["antigen_chain"].values[0].replace(" | ", "")
        if "protein" not in antigen_type:
            print('IGNORING:', antigen_type)
            continue

        antibody_chains = df_sabdab_prot[["Hchain", "Lchain"]].values[0]

        PROT_DIR = f'{EXPOSED_DIR}/{prot["idcode"]}/'
        if not os.path.isdir(PROT_DIR):
            os.mkdir(PROT_DIR)

        (
            ab_ag_complex,
            cdrchain_ag_complex,
            cdr_ag_complex,
            antibody,
            antigen,
        ) = get_systems(
            prot["idcode"],
            prot["cdr_atoms"],
            antigen_chains,
            prot["chainID"],
            antibody_chains,
            prot["cdr"],
        )
        if not antigen or not antigen.atoms:
            print('IGNORING:', prot)
            continue

        atom_res = pqr2xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr, prot["cdr"])

        run_nanoshaper(
            [ab_ag_complex, cdrchain_ag_complex, cdr_ag_complex, antibody, antigen]
        )

        (
            (ag_ab_interface, ag_ab_interface_res),
            (ag_cdrchain_interface, ag_cdrchain_interface_res),
            (ag_cdr_interface, ag_cdr_interface_res),
            (ab_ag_interface, ab_ag_interface_res),
        ) = get_interfaces(atom_res)

        prot["ag_ab_interface"] = ag_ab_interface
        prot["ag_cdrchain_interface"] = ag_cdrchain_interface
        prot["ag_cdr_interface"] = ag_cdr_interface
        prot["ab_ag_interface"] = ab_ag_interface

        prot["ag_ab_interface_res"] = ag_ab_interface_res
        prot["ag_cdrchain_interface_res"] = ag_cdrchain_interface_res
        prot["ag_cdr_interface_res"] = ag_cdr_interface_res
        prot["ab_ag_interface_res"] = ab_ag_interface_res

        new_df.append(prot.to_dict())

    print(c)

    df = pd.DataFrame(new_df)
    df.to_pickle("../data/epitope_buried.pickle")
