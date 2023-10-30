import pandas as pd
import os
import re 
import numpy as np
from  C_functs import Pdist_C,getIndex

curdir = os.path.dirname(__file__)


def get_sabdab_details():
    f = f"{curdir}/../structures/sabdab_summary_90_May2022.tsv"#f"{curdir}/../structures/sabdab_summary_90.tsv"
    df = pd.read_csv(f, sep="\t").drop_duplicates()
    df = df.query(
        "antigen_type == antigen_type and antigen_type.str.contains('protein')"
    )
    df = df.query("Hchain != Lchain")
    df = df.query("Hchain == Hchain and Lchain == Lchain")
    print(f"SabDab\nAbAg Complexes: {len(df)}\nPDB Files: {len(set(df.pdb.values))}\n")
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
                # r = str(float(line[63:70]) + 1.4) +'\n'
                newline = f"{x} {y} {z} {r}"
                xyzr.append(newline)
                coord = (x, y, z)
                atoms[anumb] = (chain, resnumb, resname, aname, cdrnumb,coord)
    with open(fout, "w") as f:
        f.write("".join(xyzr))
    return atoms


def loadTriang(annotatedTriangName):   
# print(annotatedTriangName)
    try:
        triangulationFile = open(annotatedTriangName,'r')
    except FileNotFoundError:
        print('Cannot find annotated triangulation file.\nBREAKING')
        print(annotatedTriangName)
        exit()
    triangulationFile.readline()
    triangulationFile.readline()
    triangulationFile.readline()
    infoLine=triangulationFile.readline()
    nverts=int(infoLine.split()[0])
    # print('number of vertices in the triangulation =',nverts)
    triangLines = triangulationFile.readlines()
    triangulationFile.close()
    vertLines = triangLines[:nverts]
    faceLines = triangLines[nverts:]
    


    return vertLines,faceLines

def getVolArea(NS_output):  
    '''
    IF NS error area is zero
    '''  
    # print(NS_output)
    matchV = re.search("(<<INFO>> Estimated volume )(\d*\.?\d+)",NS_output)
    matchA = re.search("(<<INFO>> Total, grid conformant, surface area is )(\d*\.?\d+)",NS_output)
    V = 0
    A = 0
    if(matchV):
        V = float(matchV.group(2))
        A = float(matchA.group(2))
    # print("AREA = ", A)

    return V,A
    


def closeVert(vert,exposedAtomCoord,thresholdDistance = 5):
    
    # not sure passing through C is interesting here
    d,_flag = Pdist_C(exposedAtomCoord,vert)
    index = np.where(d<=thresholdDistance)[0]
    if (index.size > 0):
        close = True
    else:
        close = False
    
    return close

def closeAtoms(atom1,atom2,thresholdDistance = 5):
    
    # not sure passing through C is interesting here
    d,_flag = Pdist_C(atom1,atom2)
    index = np.where(d<=thresholdDistance)[0]
    atom1_indexes,atom2_indexes = getIndex(_flag,index,atom2.shape[0])
    return atom1_indexes,atom2_indexes