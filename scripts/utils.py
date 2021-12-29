import pandas as pd
import os

curdir = os.path.dirname(__file__)


def get_sabdab_details():
    f = f"{curdir}/../structures/sabdab_summary_all.tsv"
    df = pd.read_csv(f, sep="\t")
    df = df.query("antigen_type == antigen_type and antigen_type.str.contains('protein')")
    return df


def read_pdb_line(line, pqr=False):
    aname = line[12:16].strip()
    anumb = int(line[5:11].strip())
    resname = line[17:21].strip()
    chain = line[21]
    resnumb = int(line[22:26])
    cdr_id = line[26].replace(" ", "")
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    if pqr:
        return chain, cdr_id, (resname, resnumb), (aname, anumb), (x, y, z), line[63:70]
    return chain, cdr_id, (resname, resnumb), (aname, anumb), (x, y, z)


def pqr2xyzr(fin, fout, cdrnumb):
    xyzr = []
    atoms = {}
    with open(fin) as f:
        for line in f:
            if line.startswith("ATOM "):
                chain, insertion_code, (resname, resnumb), (aname, anumb), (x, y, z) = read_pdb_line(
                    line
                )
                if insertion_code:
                    resnumb = f"{resnumb}{insertion_code}"
                r = line[63:70]
                newline = f"{x} {y} {z} {r}"
                xyzr.append(newline)
                atoms[anumb] = (chain, resnumb, resname, aname, cdrnumb)
    with open(fout, "w") as f:
        f.write("".join(xyzr))
    return atoms
