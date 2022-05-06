import os
import pandas as pd
from utils import read_pdb_line, get_sabdab_details, pqr2xyzr
import subprocess

from constants import EXPOSED_DIR, RAW_STRUCTURES_DIR


def iterate_prots(df_fname):
    # b, e = 8500, 9000
    df = pd.read_pickle(df_fname)  # .iloc[8500:]
    # print(b, e)
    print(f"CDRs: {len(df)}\tAbAg: {len(df)/6}")
    for i, row in df.iterrows():
        yield row


class System:
    def __init__(self, fname, atoms, prot_dir):
        self.fname = fname
        self.atoms = atoms
        self.indices = {}
        self.surface_atoms = []
        self.prot_dir = prot_dir
        self.fpdb = f"{prot_dir}{fname}.pdb"
        self.fpqr = f"{prot_dir}{fname}.pqr"
        self.fxyzr = f"{prot_dir}{fname}.xyzr"
        self.fexposed = f"{prot_dir}{fname}_exposed.txt"

    def trim_xyzr(self, complex_pqr, complex_xyzr):
        self_xyzr = []
        idx = 0
        with open(complex_pqr) as f_pqr, open(complex_xyzr) as f_xyzr:
            for line_pqr, line_xyzr in zip(f_pqr, f_xyzr):
                chain, (_, resnumb), (_, anumb), (x, y, z), radius = read_pdb_line(
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
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        return proc.stdout.decode("ascii")
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"-> {command} <- did not run successfully\nMessage: {e.stdout.decode('ascii')}"
        )


def get_systems(idcode, cdr_atoms, antigen_chains, cdr_chain, ab_chains, cdr_numb):
    complex_pdb = []
    ab_cdr_atoms = []
    chain_res = {}
    print(f"{RAW_STRUCTURES_DIR}/{idcode}.pdb")
    with open(f"{RAW_STRUCTURES_DIR}/{idcode}.pdb") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            (
                chain,
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

    f_complex = f"{idcode}_complex_{chains}_{antigen_chains}"  # antibody + antigen
    f_chaincomplex = f"{idcode}_chaincomplex_{cdr_chain}_{antigen_chains}"  # one antibody chain + antigen
    f_cdrcomplex = f"{idcode}_cdrcomplex_{cdr_chain}{cdr_numb}_{antigen_chains}"  # one CDR + antigen
    f_antibody = f"{idcode}_antibody_{chains}"  # only antibody
    f_antigen = f"{idcode}_antigen_{antigen_chains}"  # only antigen

    prot_dir = f"{EXPOSED_DIR}/{idcode}/"

    ab_ag_complex = System(
        f_complex,
        {
            chain: list(set(chain_res[chain]))
            for chain in antigen_chains + cdr_chain + otherab_chain
        },
        prot_dir,
    )
    cdrchain_ag_complex = System(
        f_chaincomplex,
        {chain: list(set(chain_res[chain])) for chain in antigen_chains + cdr_chain},
        prot_dir,
    )
    cdr_ag_complex = System(
        f_cdrcomplex,
        dict(
            {chain: list(set(chain_res[chain])) for chain in antigen_chains},
            **{cdr_chain: ab_cdr_atoms},
        ),
        prot_dir,
    )
    # cdrchain = System(f_cdrchain, {cdr_chain: ab_cdrchain_res})
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

    pdb2pqr_err = False
    if not os.path.isfile(ab_ag_complex.fpqr):
        with open(ab_ag_complex.fpdb, "w") as f:
            f.write("".join(complex_pdb))

        cmd = f"""
        pdb2pqr30 --ff=AMBER --keep-chain --noopt --nodebump {ab_ag_complex.fpdb} {ab_ag_complex.fpqr};
        if [ -f {ab_ag_complex.fpqr} ];
        then 
            grep "ATOM" {ab_ag_complex.fpqr} > {prot_dir}/tmp.pqr;
            mv {prot_dir}/tmp.pqr {ab_ag_complex.fpqr};
        fi
        """
        stdout = run_command(cmd)

        # check if pdb2pqr crashed
        if not os.path.isfile(ab_ag_complex.fpqr):
            pdb2pqr_err = stdout
    else:
        print(f"SKIPPING pdb2pqr run of {ab_ag_complex.fpdb}")

    return (
        ab_ag_complex,
        cdrchain_ag_complex,
        cdr_ag_complex,
        antibody,
        antigen,
        pdb2pqr_err,
    )


def run_nanoshaper(systems: list, recalc: bool):
    for system in systems:
        print(system.fname)
        system.trim_xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr)

        if os.path.isfile(system.fexposed) and not recalc:
            print("SKIPPED", system.fexposed)
            system.get_surface_atoms()
            continue

        cmd = f"""sed "s/__XYZRFILENAME__/{system.fxyzr.split('/')[-1]}/" template.prm > {system.prot_dir}conf.prm"""
        os.system(cmd)

        print("running NanoShaper")
        cmd = f"docker run -i --cpus=8 --mount type=bind,source={system.prot_dir},target=/App registry-gitlab.iit.it/sdecherchi/nanoshaper:0.7.8 conf.prm"
        run_command(cmd)

        os.system(
            f"cp {system.prot_dir}exposedIndices.txt {system.fexposed}; rm -f {system.prot_dir}exposedIndices.txt"
        )
        system.get_surface_atoms()
    return


def get_interfaces(atom_res):
    ag_ab_interface = antigen.surface_atoms - ab_ag_complex.surface_atoms
    ag_cdrchain_interface = antigen.surface_atoms - cdrchain_ag_complex.surface_atoms
    ag_cdr_interface = antigen.surface_atoms - cdr_ag_complex.surface_atoms
    ab_ag_interface = antibody.surface_atoms - ab_ag_complex.surface_atoms
    ag_interface = antigen.surface_atoms

    ag_ab_interface_res = [atom_res[i] for i in ag_ab_interface]
    ag_cdrchain_interface_res = [atom_res[i] for i in ag_cdrchain_interface]
    ag_cdr_interface_res = [atom_res[i] for i in ag_cdr_interface]
    ab_ag_interface_res = [atom_res[i] for i in ab_ag_interface]
    ag_interface_res = [atom_res[i] for i in ag_interface]

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
        (ag_interface, ag_interface_res),
    )


if __name__ == "__main__":
    df_sabdab_all = get_sabdab_details()

    f_epitopes = "../data/cdr_epitope.pickle"

    new_df = []
    c = 0
    pdb2pqr_fails = "PDB2PQR_FAILS"
    os.system(f"rm -f {pdb2pqr_fails}")
    for prot in iterate_prots(f_epitopes):
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
            print("IGNORING:", antigen_type)
            continue

        antibody_chains = df_sabdab_prot[["Hchain", "Lchain"]].values[0]

        if prot["idcode"] == "4cc8" and antigen_chains == "G":
            print("IGNORING 4cc8")
            continue

        prot_dir = f'{EXPOSED_DIR}/{prot["idcode"]}/'
        if not os.path.isdir(prot_dir):
            os.mkdir(prot_dir)

        (
            ab_ag_complex,
            cdrchain_ag_complex,
            cdr_ag_complex,
            antibody,
            antigen,
            pdb2pqr_err,
        ) = get_systems(
            prot["idcode"],
            prot["cdr_atoms"],
            antigen_chains,
            prot["chainID"],
            antibody_chains,
            prot["cdr"],
        )

        if pdb2pqr_err:
            print(pdb2pqr_err)
            with open(pdb2pqr_fails, "a") as f:
                f.write(f"{prot['idcode']}\n")
            continue

        if not antigen.atoms:
            print(antigen)
            print(antigen.atoms)
            print("IGNORING:", prot)
            exit()
            continue

        atom_res, recalc = pqr2xyzr(
            ab_ag_complex.fpqr, ab_ag_complex.fxyzr, prot["cdr"]
        )

        run_nanoshaper(
            [ab_ag_complex, cdrchain_ag_complex, cdr_ag_complex, antibody, antigen],
            recalc,
        )

        (
            (ag_ab_interface, ag_ab_interface_res),
            (ag_cdrchain_interface, ag_cdrchain_interface_res),
            (ag_cdr_interface, ag_cdr_interface_res),
            (ab_ag_interface, ab_ag_interface_res),
            (ag_interface, ag_interface_res),
        ) = get_interfaces(atom_res)

        prot["ag_ab_interface"] = ag_ab_interface  # Ag residues interfacing with the Ab
        prot[
            "ag_cdrchain_interface"
        ] = ag_cdrchain_interface  # Ab residues interfacing with one Ab chain
        prot[
            "ag_cdr_interface"
        ] = ag_cdr_interface  # Ag residues interfacing with one CDR
        prot["ab_ag_interface"] = ab_ag_interface  # Ab residues interfacing with the Ag
        prot["ag_interface"] = ag_interface

        prot["ag_ab_interface_res"] = ag_ab_interface_res
        prot["ag_cdrchain_interface_res"] = ag_cdrchain_interface_res
        prot["ag_cdr_interface_res"] = ag_cdr_interface_res
        prot["ab_ag_interface_res"] = ab_ag_interface_res
        prot["ag_interface_res"] = ag_interface_res

        new_df.append(prot.to_dict())

    print(c)

    df = pd.DataFrame(new_df)
    df.to_pickle("../data/epitope_buried_.pickle")
