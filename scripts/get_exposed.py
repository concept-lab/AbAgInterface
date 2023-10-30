#TODO 
# Run saving all possibly necessary info: 
# 1. anumb interface for each structure
# 2. full info including coordinaes of each atom
# 3. Triangulation with vertices for each system

import os
import pandas as pd
from utils import getVolArea, loadTriang, read_pdb_line, get_sabdab_details, pqr2xyzr,closeVert,closeAtoms
import subprocess
import numpy as np
from time import time
from datetime import datetime


colorMap = True
thresholD = 4
useInterior = True

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
EXPOSED_DIR = f"{CUR_DIR}/../structures/exposed/"
print("EXPOSED DIR= ",EXPOSED_DIR)
WORKDIR = f'{CUR_DIR}/../tempNS_prova/'
print("WORKING DIR= ",WORKDIR)
print("/n Build colored triangulation= ",colorMap)

print("Post filtering threshold between vertices and exposed atoms= ",thresholD)
print("Extend to vertices of close atoms= ",useInterior)

input("\nCONTINUE?")


def iterate_prots(df):
    # df = pd.read_pickle(df_fname)  # .iloc[::-1]
    for i, row in df.iterrows():
        yield row
    # exit()


class System:
    def __init__(self, fname, atoms,type,chainName=None):
        self.fname = fname
        self.atoms = atoms
        
        self.indices = {} #map between index in the xyzr file (NS indicization) and atom number in PQR file, used as absolute reference and to connect to res
        self.surface_atoms = []
        self.fpdb = f"{PROT_DIR}{fname}.pdb"
        self.fpqr = f"{PROT_DIR}{fname}.pqr"
        self.fxyzr = f"{PROT_DIR}{fname}.xyzr"
        self.fexposed = f"{PROT_DIR}{fname}_exposed.txt"
        self.fexposed9 = f"{PROT_DIR}{fname}_exposed9.txt"
        self.fexposed100 = f"{PROT_DIR}{fname}_exposed100.txt"
        self.ftriang = f"{PROT_DIR}{fname}_annotated.off"
        self.fout = f"{PROT_DIR}{fname}_NS.out"

        self.vertLines = [] 
        self.faceLines = [] 
        self.type = type
        self.area = None
        # self.surface_atoms_all = set()
        self.allAtoms = set()

        self.chainName = chainName
        self.surface_atoms9 = []
        self.surface_atoms100 = []

        # self.xyzr=[]
        

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
                    # self.xyzr.append(line_xyzr)
                    self.indices[idx] = anumb
                    self.allAtoms.add(anumb)
                    idx += 1

                
        if (self.type=='complex'):
            cmd = f'cp {self.fxyzr} {WORKDIR}/NS_input.xyzr'
            run_command(cmd)
            # with open(self.fxyzr, "w") as f:
            #     f.write("".join(self_xyzr)) #permanently store only complex xyzr which is needed
        else:
            with open(f'{WORKDIR}/NS_input.xyzr', "w") as f:
                f.write("".join(self_xyzr))

    def get_surface_atoms(self,R=1.4):
        if (R==9):
            try:
                with open(self.fexposed9) as f:
                    for line in f:
                        idx = line.strip()
                        if idx:
                            idx = int(idx)
                            try:
                                anumb = self.indices[idx]
                            except:
                                print("KEY ERROR in EXPOSED ATOMS.. Probably due to wrong exposed file")
                                print(self.fexposed)
                                print(self.fout)
                                input('continue?')
                                self.surface_atoms9 = set()
                                break
                            self.surface_atoms9.append(anumb)
                self.surface_atoms9 = set(self.surface_atoms9)
            except:
                print("could not fetch exposed at rp=9")
                self.surface_atoms9=set()
        elif (R==100):
            try:
                with open(self.fexposed100) as f:
                    for line in f:
                        idx = line.strip()
                        if idx:
                            idx = int(idx)
                            try:
                                anumb = self.indices[idx]
                            except:
                                print("KEY ERROR in EXPOSED ATOMS.. Probably due to wrong exposed file")
                                print(self.fexposed)
                                print(self.fout)
                                input('continue?')
                                self.surface_atoms100 = set()
                                break
                            self.surface_atoms100.append(anumb)
                self.surface_atoms100 = set(self.surface_atoms100)
            except:
                print("could not fetch exposed at rp=100")
                self.surface_atoms100=set()
        else:
            with open(self.fexposed) as f:
                for line in f:
                    idx = line.strip()
                    if idx:
                        idx = int(idx)
                        try:
                            anumb = self.indices[idx]
                        except:
                            print("KEY ERROR in EXPOSED ATOMS.. Probably due to wrong exposed file")
                            print(self.fexposed)
                            exit()
                        #     print(self.fout)
                            # input('continue?')
                            # self.surface_atoms = []
                            # break
                        self.surface_atoms.append(anumb)
            self.surface_atoms = set(self.surface_atoms)

        # vertLines = loadTriang(f'{WORKDIR}triangulatedSurf.off')
        # if(self.type=='antigen' or self.type=='antibody'):
        #     self.vertLines=vertLines
        
        # for line in vertLines:
        #     idx = line.split()[6]
        #     # if idx:
        #     idx = int(idx)
        #     anumb = self.indices[idx]
        #     self.surface_atoms_all.add(anumb)


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
    # print(f"{base}/{idcode}.pdb")
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
    f_chain= f"{idcode}_chain_{cdr_chain}"
    f_cdrcomplex = f"{idcode}_cdrcomplex_{cdr_chain}{cdr_numb}_{antigen_chains}"
    f_cdr= f"{idcode}_cdr_{cdr_chain}{cdr_numb}"
    f_antibody = f"{idcode}_antibody_{chains}"
    f_antigen = f"{idcode}_antigen_{antigen_chains}"
    
    
    ab_ag_complex = System(
        f_complex,
        {
            chain: list(set(chain_res[chain]))
            for chain in antigen_chains + cdr_chain + otherab_chain
        }, type="complex"
    )
    cdrchain_ag_complex = System(
        f_chaincomplex,
        {chain: list(set(chain_res[chain])) for chain in antigen_chains + cdr_chain},
        type="cdrchainComplex"
    )

    cdrchain = System(
        f_chain,
        {chain: list(set(chain_res[chain])) for chain in cdr_chain},
        type="cdrchain"
    )

    cdr_ag_complex = System(
        f_cdrcomplex,
        dict(
            {chain: list(set(chain_res[chain])) for chain in antigen_chains},
            **{cdr_chain: ab_cdr_atoms},
        ),type="cdrComplex"
    )

    cdr = System(
        f_cdr,
        {cdr_chain: ab_cdr_atoms}
        ,type="cdr"
    )
    # cdrchain = System(f_cdrchain, {cdr_chain: ab_cdrchain_res})
    antibody = System(
        f_antibody,
        {chain: list(set(chain_res[chain])) for chain in cdr_chain + otherab_chain},
        type="antibody"
    )
    antigen = System(
        f_antigen, {chain: list(set(chain_res[chain])) for chain in antigen_chains},
        type="antigen",
        chainName = f"antigen_{antigen_chains}"
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
            # print("pdb2pqr FAILED for "+ab_ag_complex.fname)
            return [None] * 7
    else:
        # print(f"SKIPPING pdb2pqr run of {ab_ag_complex.fpdb}")
        pass

    return (ab_ag_complex, cdrchain_ag_complex, cdr_ag_complex, antibody, antigen,cdrchain,cdr)


def run_nanoshaper(systems: list):
    for system in systems:
        # print(system.fname)
        system.trim_xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr) #stores internally indexes of atoms of each (sub)system and copies to NS temp folder NS input
        
        if os.path.isfile(system.fout):
            fNSout = open(system.fout)
            out = fNSout.readlines()
            _V,A = getVolArea(str(out))
            system.area = A
            if(os.path.isfile(system.fexposed) and system.area>0):   
                #it could be that there are some previously bad exposed files, which were wronlgy copied from previous run even after NS failure..
                # now mv and not copy to prevent this scenario    
                system.get_surface_atoms() #EXPOSED AT NORMAL RADIUS
                
            else:
                #In interface computations if surface atoms has not been filled it should raise error (list operation "-" does not exists)
                print("AREA= ",system.area )
                print("NS error:", system.fout)
            if(system.type == "antigen"):# or system.type == "complex"):
                ok=0
                # print('here')
                if(os.path.isfile(system.fexposed9)):
                    system.get_surface_atoms(R=9)
                    ok+=1
                if(os.path.isfile(system.fexposed100)):
                    system.get_surface_atoms(R=100)
                    ok+=1
                if ok<2:
                    
                    # print('NS new runs for larger rp: ')
                    # print(system.fname)
                    try:
                        out = subprocess.check_output(['./NanoShaper', 'template9.prm'],cwd=WORKDIR)
                        os.system(f"mv {WORKDIR}exposedIndices.txt {system.fexposed9}")
                    except subprocess.CalledProcessError as grepexc:
                        print("CANNOT GET EXPOSED AT PROBE RADIUS=9")                                                                                                   
                        print ("error code", grepexc.returncode, grepexc.output)
                    system.get_surface_atoms(R=9)
                    try:
                        out = subprocess.check_output(['./NanoShaper', 'template100.prm'],cwd=WORKDIR)
                        os.system(f"mv {WORKDIR}exposedIndices.txt {system.fexposed100}")
                    except subprocess.CalledProcessError as grepexc:    
                        print("CANNOT GET EXPOSED AT PROBE RADIUS=100")                                                                                                 
                        print ("error code", grepexc.returncode, grepexc.output)
                    system.get_surface_atoms(R=100)
                else:
                    # print(system.fname)
                    # print("rp =9 and rp =100 already computed")
                    pass
                
                
            
            continue

        # print('NS new runs: ')
        # print(system.fname)
        
            
        NS_error=False
        try:
            out = subprocess.check_output(['./NanoShaper', 'template.prm'],cwd=WORKDIR) #WORKDIR is where all NS runs happen, the output is copied in PROT_DIR 
        except subprocess.CalledProcessError as grepexc:                                                                                                   
            print ("error code", grepexc.returncode, grepexc.output)
            out = grepexc.output
            NS_error = True
        if(not NS_error):
            #Happens onlyfor cdrs in case which is not important for exposed indices whch are taken from the complex
            os.system(
                f"mv {WORKDIR}exposedIndices.txt {system.fexposed}"
            )

            # print(out)
            system.get_surface_atoms()
            # print(system.surface_atoms)

        else:
            print("NS error", system.fname)
            pass

        # print('here')
        ############ STORE OUTPUT ########
        fout = open(system.fout,'w')
        fout.write(str(out))
        fout.close()
        _V,A = getVolArea(str(out))
        system.area = A
        ######### SAVE TRIANGULATIONS with normals snf closer atom index annotation only of antigen and antib
        if(system.type == "antigen" or system.type == "antibody"):
            cmd = f'mv {WORKDIR}triangulatedSurf.off {system.ftriang}'
            run_command(cmd)
        else:
            pass

        #########

        if(system.type == "antigen"):# or system.type == "complex"):
            ok=0
            if(os.path.isfile(system.fexposed9)):
                system.get_surface_atoms(R=9)
                ok+=1
            if(os.path.isfile(system.fexposed100)):
                system.get_surface_atoms(R=100)
                ok+=1
            if ok<2:
                
                # print('NS new runs for larger rp: ')
                # print(system.fname)
                try:
                    out = subprocess.check_output(['./NanoShaper', 'template9.prm'],cwd=WORKDIR)
                    os.system(f"mv {WORKDIR}exposedIndices.txt {system.fexposed9}")
                except subprocess.CalledProcessError as grepexc:
                    print("CANNOT GET EXPOSED AT PROBE RADIUS=9")                                                                                                   
                    print ("error code", grepexc.returncode, grepexc.output)
                system.get_surface_atoms(R=9)
                try:
                    out = subprocess.check_output(['./NanoShaper', 'template100.prm'],cwd=WORKDIR)
                    os.system(f"mv {WORKDIR}exposedIndices.txt {system.fexposed100}")
                except subprocess.CalledProcessError as grepexc:    
                    print("CANNOT GET EXPOSED AT PROBE RADIUS=100")                                                                                                 
                    print ("error code", grepexc.returncode, grepexc.output)
                system.get_surface_atoms(R=100)
            else:
                # print("already computed")
                pass

        

    return



def get_interfaces(atom_res):
    '''IMPORTANT: setting by hand to 0 areas if no exposed residues (actually there could still be cotact area even if no exposed interface atoms..)'''

    #TODO implementare contatore interno che distingue che sto nella stessa struttura ma sto cambiando cdr. --> fare vertex map distinta per ogni cdr
    # Credo no problem per la catena visto che genera tutta una nuova triangolazione dell'antigene

    noColor = [128, 128, 128] #none
    violet = [255,0,255] #all anitgen antibody interface VIOLET
# palette_extended[1,:] = [0,204,255] #all anitgen antibody interface BLUE
    yellow = [255,204,0] #CDR chain YELLOW
    red = [255,0,0] #CDR  RED


    failure=None
    gotFailure = False

    # Comment: is a set substraction
    ag_ab_interface = antigen.surface_atoms - ab_ag_complex.surface_atoms
    ag_cdrchain_interface = antigen.surface_atoms - cdrchain_ag_complex.surface_atoms #CONTACT FROM THE ANTIGEN SIDE
    ag_cdr_interface = antigen.surface_atoms - cdr_ag_complex.surface_atoms
    ab_ag_interface = antibody.surface_atoms - ab_ag_complex.surface_atoms

    ag_interface = antigen.surface_atoms

    ag_ab_interface_res = [atom_res[i] for i in ag_ab_interface]
    ag_cdrchain_interface_res = [atom_res[i] for i in ag_cdrchain_interface]
    ag_cdr_interface_res = [atom_res[i] for i in ag_cdr_interface]
    ab_ag_interface_res = [atom_res[i] for i in ab_ag_interface]

    ag_interface_res = [atom_res[i] for i in ag_interface]

    ag_interface9 = antigen.surface_atoms9
    ag_interface100 = antigen.surface_atoms100
    antigen_res9 = [atom_res[i] for i in ag_interface9]
    antigen_res100 = [atom_res[i] for i in ag_interface100]
    # ag_ab_interface9 = antigen.surface_atoms9 - ab_ag_complex.surface_atoms9
    # ag_ab_interface100 = antigen.surface_atoms100 - ab_ag_complex.surface_atoms100
    # ag_ab_interface_res9 = [atom_res[i] for i in ag_ab_interface9]
    # ag_ab_interface_res100 = [atom_res[i] for i in ag_ab_interface100]

    # print(set(ag_ab_interface_res))
    # print(set(ag_cdrchain_interface_res))
    # print(set(ag_cdr_interface_res))
    # print(set(ab_ag_interface_res))
    # exit()


    ################# NEW ################
    # NOW GET ALSO AREA OF CONTACT 

    # following to avoid to run NS and use out file stores 
    if(ab_ag_complex.area is None):
        fout= open(ab_ag_complex.fout,'r')
        _V,A = getVolArea(fout.readline())
        ab_ag_complex.area = A
    if(antigen.area is None):
        fout= open(antigen.fout,'r')
        _V,A = getVolArea(fout.readline())
        antigen.area = A
    if(antibody.area is None):
        fout= open(antibody.fout,'r')
        _V,A = getVolArea(fout.readline())
        antibody.area = A
    if(cdrchain_ag_complex.area is None):
        fout= open(cdrchain_ag_complex.fout,'r')
        _V,A = getVolArea(fout.readline())
        cdrchain_ag_complex.area = A
    if(cdr_ag_complex.area is None):
        fout= open(cdr_ag_complex.fout,'r')
        _V,A = getVolArea(fout.readline())
        cdr_ag_complex.area = A

    ag_ab_interface_AREA = 0.5*abs(antigen.area - (ab_ag_complex.area - antibody.area)) #--> acutally is an average, assuming that area both interfaces equal 
    # ab_ag_interface_AREA = 0.5*abs(antibody.area - (ab_ag_complex.area - antigen.area)) #IDENTICAL ONLY if no cavities or cavities not filled
                                                                    # area of cdrchain not complex
    # print((ab_ag_complex.area-cdrchain_ag_complex.area ))
    # print(antibody.area)
    # ag_cdrchain_AREA = antigen.area - (cdrchain_ag_complex.area  - (antibody.area - (ab_ag_complex.area-cdrchain_ag_complex.area )))
    ag_cdrchain_AREA = 0.5*abs(antigen.area - (cdrchain_ag_complex.area - cdrchain.area))
    # ag_cdr_AREA = antigen.area - (cdr_ag_complex.area  - (antibody.area - (ab_ag_complex.area-cdr_ag_complex.area )))
    ag_cdr_AREA= 0.5*abs(antigen.area - (cdr_ag_complex.area - cdr.area)) 

    # print("antigenAREA=",antigen.area)
    # print("cdrchain complex area=",cdrchain_ag_complex.area)
    # print("cdrchain area=",cdrchain.area)
    # print("cdr complex area=",cdr_ag_complex.area)
    # print("cdr area=",cdr.area)
    # print("ag-cdr area=",ag_cdr_AREA)
    # input()

    # SELECT VERTICES CLOSE TO CONTACT AREA FOR EACH SYSTEM 
    # print(antigen.ftriang)
    antigen.vertLines,antigen.faceLines = loadTriang(antigen.ftriang)
    # print(antibody.ftriang)
    antibody.vertLines,antibody.faceLines = loadTriang(antibody.ftriang) 
    
    

    atomResMap = {}
    for anumb,r in atom_res.items():
        if((r[0],r[1]) in atomResMap):
            atomResMap[(r[0],r[1])].add(anumb)#(chain,ID)-->atom list
        else:
            atomResMap[(r[0],r[1])]=set([anumb])

    
    if(ag_ab_interface_res):

        ag_ab_interfaceAtoms_coord = np.array([r[5] for r in ag_ab_interface_res]) #not res but full coordinates.. the name confounds a lot
        ab_ag_interfaceAtoms_coord = np.array([r[5] for r in ab_ag_interface_res])
        
        ag_inside_atoms = antigen.allAtoms-ag_ab_interface #rather than inside these are all the others not considered exposed at the interface
        ab_inside_atoms = antibody.allAtoms-ab_ag_interface 
    
        ag_inside_atoms = np.array(list(ag_inside_atoms))
        ab_inside_atoms = np.array(list(ab_inside_atoms))
        ag_inside_atoms_coord = np.array([atom_res[i][5] for i in ag_inside_atoms])
        ab_inside_atoms_coord = np.array([atom_res[i][5] for i in ab_inside_atoms])
        
        indx_ag,indx_interface = closeAtoms(ag_inside_atoms_coord,ag_ab_interfaceAtoms_coord,thresholdDistance=5)
        indx_ab,indx_interface = closeAtoms(ab_inside_atoms_coord,ab_ag_interfaceAtoms_coord,thresholdDistance=5)
        
        insideInterface_ag = set(ag_inside_atoms[indx_ag])
        insideInterface_ab = set(ab_inside_atoms[indx_ab])
        if useInterior:
            mapRES_ag_ab_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ag_ab_interface_res]) | insideInterface_ag #first set is list of atom numbers res based, second set is atom based
        # print(mapRES_ag_ab_interface)
        # input()
         # mapRES_ab_ag_interface = ab_ag_interface | insideInterface_ab
        # mapRES_ag_ab_interface = ag_ab_interface | insideInterface_ag
            mapRES_ab_ag_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ab_ag_interface_res]) | insideInterface_ab
        else:
            mapRES_ag_ab_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ag_ab_interface_res])
            mapRES_ab_ag_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ab_ag_interface_res])
       
    else:
        mapRES_ag_ab_interface = set()
        mapRES_ab_ag_interface = set()
        print('NO AB AG CONTACT!'+ab_ag_complex.fname)
        failure='NO AB AG CONTACT!'+ab_ag_complex.fname
        gotFailure = True
        ag_ab_interface_AREA = 0
        ab_ag_interface_AREA = 0
        ag_cdrchain_AREA = 0 
        ag_cdr_AREA = 0

    noContact = False
    if(ag_cdrchain_interface_res):
        mapRES_cdrchain_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ag_cdrchain_interface_res]) #| insideInterface_ag
    else:
        mapRES_cdrchain_interface = set()
        print("THIS CDRCHAIN IS NOT IN CONTACT :"+cdrchain_ag_complex.fname)
        ag_cdrchain_AREA = 0 
        noContact = True
        if(not gotFailure):
            failure='THIS CDRCHAIN IS NOT IN CONTACT: '+cdrchain_ag_complex.fname
            gotFailure = True
            noContact = True
    # mapRES_cdrchain_interface = ag_cdrchain_interface #| insideInterface_ag
    if(ag_cdr_interface_res):
        mapRES_cdr_interface = set.union(*[atomResMap[(r[0],r[1])] for r in ag_cdr_interface_res]) #| insideInterface_ag
    else:
        mapRES_cdr_interface = set()
        print("THIS CDR IS NOT IN CONTACT: "+cdr_ag_complex.fname)
        ag_cdr_AREA = 0
        if(not gotFailure):
            failure="THIS CDR IS NOT IN CONTACT: "+ cdr_ag_complex.fname
        noContact = True
    # mapRES_cdr_interface =ag_cdr_interface #| insideInterface_ag

    
    
    
    vertLine=set()
    # printVert=[]
    for ind,line in enumerate(antigen.vertLines):
        l=line.split()[0:3]
        l.append('\n')
        # print(l)
        # input()
        # printVert.append(' '.join(l))
        indAtomClose = int(line.split()[6]) #referred to index in antigen.xyzr
        anumb = antigen.indices[indAtomClose]
        x = float(line.split()[0])
        y = float(line.split()[1])
        z = float(line.split()[2])
        coordVert = np.array([[x,y,z]])
        if anumb in mapRES_ag_ab_interface:
            #INVERSION TO MAKE THEM CLOSER TO OTHER CONTACT (INDEX OF ATOM CLOSE TO THE OPPOSITE SURFACE)
            isClose = closeVert(coordVert,ab_ag_interfaceAtoms_coord,thresholdDistance=thresholD)
            if (isClose):
                vertLine.add(ind)
    
    # print(vertLine)
    #contact patch on antigen
    if(len(vertLine)>0):
        if os.path.isfile(f'{PROT_DIR}'+antibody.fname+'_'+antigen.chainName+'_AGcontactPatch.off'):
            pass
        else:
            contactPatchFile = open(f'{PROT_DIR}'+antibody.fname+'_'+antigen.chainName+'_AGcontactPatch.off','w')
            
            contacts=[]
            for line in antigen.faceLines:
                faceConnections = [int (i) for i in line.split()[1:]]
                # print(faceConnections)
                # input()
                if (len(set(faceConnections).intersection(vertLine))>0): # 2 or more vertices are tagged
                    contacts.append(line)
            
            contactPatchFile.write("OFF\n")
            contactPatchFile.write("%d\t%d\t0\n"%(len(antigen.vertLines),len(contacts)))
            # contactPatchFile.writelines(printVert) #other columns contain normals and closest atom infos
            contactPatchFile.writelines(" ".join([*s.split()[0:3],'\n']) for s in antigen.vertLines)
            contactPatchFile.writelines(contacts)
            contactPatchFile.close()
            
            # vert map for personalised coloring
        if colorMap and (noContact==False): #avoid creating unnecessary files
            # vertFile = open(f'{PROT_DIR}'+cdr_ag_complex.fname+'_vertMapResBased.txt','w')
            coloredTriang = open(f'{PROT_DIR}'+cdr_ag_complex.fname+'_coloredTriang.off','w')
            coloredTriang.write("# Red CDR, Yellow whole chain (H or L), Violet remaining contact\n")
            coloredTriang.write("COFF\n")
            coloredTriang.write("%d\t%d\t0\n"%(len(antigen.vertLines),len(antigen.faceLines)))
            for ind,line in enumerate(antigen.vertLines):
                flag=0
                color = noColor
                if ind in vertLine:
                    flag=1
                    color = violet
                    indAtomClose = int(antigen.vertLines[ind].split()[6]) #referred to index in antigen.xyzr
                    anumb = antigen.indices[indAtomClose]
                    if anumb in mapRES_cdrchain_interface: #these are by construction subsets of the above
                        flag = 2 
                        color = yellow
                        if anumb in mapRES_cdr_interface: # same
                            flag = 3
                            color = red
                coloredTriang.write(" ".join(line.split()[0:3]))
                coloredTriang.write(f" {int(color[0])} ")
                coloredTriang.write(f"{int(color[1])} ") 
                coloredTriang.write(f"{int(color[2])} ") 
                coloredTriang.write(f"255\n")
            coloredTriang.writelines(antigen.faceLines)
            coloredTriang.close()
            #     vertFile.write(str(flag)+'\n')
            # vertFile.close()
    
    # for line in antigen.vertLines:
    #     indAtomClose = int(line.split()[6]) #referred to index in antigen.xyzr
    #     anumb = antigen.indices[indAtomClose]
    #     x = float(line.split()[0])
    #     y = float(line.split()[1])
    #     z = float(line.split()[2])
    #     coordVert = np.array([[x,y,z]])
    #     flag = 0
    #     if anumb in mapRES_ag_ab_interface:
    #         isClose = closeVert(coordVert,ag_ab_interfaceAtoms_coord,thresholdDistance=5)
    #         if (isClose):
    #             flag = 1
    #             if anumb in mapRES_cdrchain_interface: #these are by construction subsets of the above
    #                 flag = 2 
    #             if anumb in mapRES_cdr_interface: # same
    #                 flag = 3
    #     vertFile.write(str(flag)+'\n')
    # vertFile.close()






    
    if os.path.isfile(f'{PROT_DIR}'+antibody.fname+'_vertMapResBased.txt'):
        pass
    else:
        vertLine=set()
        # printVert=[]
        for ind,line in enumerate(antibody.vertLines):
            l=line.split()[0:3]
            l.append('\n')
            # print(l)
            # input()
            # printVert.append(' '.join(l))
            indAtomClose = int(line.split()[6]) #referred to index in antigen.xyzr
            anumb = antibody.indices[indAtomClose]
            x = float(line.split()[0])
            y = float(line.split()[1])
            z = float(line.split()[2])
            coordVert = np.array([[x,y,z]])
            if anumb in mapRES_ab_ag_interface:
                #INVERSION TO MAKE THEM CLOSER TO OTHER CONTACT (INDEX OF ATOM CLOSE TO THE OPPOSITE SURFACE)
                isClose = closeVert(coordVert,ag_ab_interfaceAtoms_coord,thresholdDistance=thresholD)
                if (isClose):
                    vertLine.add(ind)

#contact patch on antibody
        if(len(vertLine)>0):
            contactPatchFile = open(f'{PROT_DIR}'+antibody.fname+'_'+antigen.chainName+'_ABcontactPatch.off','w')
            contacts=[]
            for line in antibody.faceLines:
                faceConnections = [int (i) for i in line.split()[1:]]
                # print(faceConnections)
                # input()
                if (len(set(faceConnections).intersection(vertLine))>0): # 2 or more vertices are tagged
                    contacts.append(line)
            
            contactPatchFile.write("COFF\n")
            contactPatchFile.write("%d\t%d\t0\n"%(len(antibody.vertLines),len(contacts)))
            # contactPatchFile.writelines(printVert) #other columns contain normals and closest atom infos
            contactPatchFile.writelines(" ".join([*s.split()[0:3],'\n']) for s in antibody.vertLines)
            contactPatchFile.writelines(contacts)
            contactPatchFile.close()
            

            if colorMap:
                # vertFile = open(f'{PROT_DIR}'+antibody.fname+'_vertMapResBased.txt','w')
                coloredTriang = open(f'{PROT_DIR}'+antibody.fname+'_coloredTriang.off','w')
                coloredTriang.write("# Red CDR, Yellow whole chain (H or L), Violet remaining contact\n")
                coloredTriang.write("COFF\n")
                coloredTriang.write("%d\t%d\t0\n"%(len(antibody.vertLines),len(antibody.faceLines)))


                for ind,line in enumerate(antibody.vertLines):
                    flag=0
                    color = noColor
                    if ind in vertLine:
                        flag=1
                        color = violet
                        # indAtomClose = int(antibody.vertLines[ind].split()[6]) #referred to index in antigen.xyzr
                        # anumb = antibody.indices[indAtomClose]
                    coloredTriang.write(" ".join(line.split()[0:3]))
                    coloredTriang.write(f" {int(color[0])} ")
                    coloredTriang.write(f"{int(color[1])} ") 
                    coloredTriang.write(f"{int(color[2])} ") 
                    coloredTriang.write(f"255\n")
                coloredTriang.writelines(antibody.faceLines)
                coloredTriang.close()
                #     vertFile.write(str(flag)+'\n')
                # vertFile.close()
        
        # for line in antibody.vertLines:
        #     indAtomClose = int(line.split()[6])
        #     anumb = antibody.indices[indAtomClose]
        #     flag = 0
        #     x = float(line.split()[0])
        #     y = float(line.split()[1])
        #     z = float(line.split()[2])
        #     coordVert = np.array([[x,y,z]])
            
        #     if anumb in mapRES_ab_ag_interface:
        #         isClose = closeVert(coordVert,ab_ag_interfaceAtoms_coord,thresholdDistance=5)
        #         if isClose:
        #             flag = 1
                    
        #     vertFile.write(str(flag)+'\n')
        

    #PROVE (on AG SIDE ONLY):
    # interfaceRES = set( [(r[0],r[1]) for r in ag_ab_interface_res] )
    # # interfaceRES9 = set( [(r[0],r[1]) for r in ag_ab_interface_res9] )
    # # interfaceRES100 = set( [(r[0],r[1]) for r in ag_ab_interface_res100] )
    # antigenRES9 = set( [(r[0],r[1]) for r in antigen_res9] )
    # antigenRES100 = set( [(r[0],r[1]) for r in antigen_res100] )
    
    # print("Intersection rp=1.4 with rp=9",interfaceRES.intersection(antigenRES9))
    # # print("Intersection rp=1.4 with rp=9 USING INTERFACES",interfaceRES.intersection(interfaceRES9))
    # print("Intersection rp=1.4 with rp=100",interfaceRES.intersection(antigenRES100))
    # print("Intersection rp=1.4 with rp=100 USING INTERFACES",interfaceRES.intersection(interfaceRES100))
    # print(interfaceRES)

    return (failure,
        (ag_ab_interface, ag_ab_interface_res,ag_ab_interface_AREA),
        (ag_cdrchain_interface, ag_cdrchain_interface_res,ag_cdrchain_AREA),
        (ag_cdr_interface, ag_cdr_interface_res,ag_cdr_AREA),
        (ab_ag_interface, ab_ag_interface_res),
        (ag_interface, ag_interface_res,antigen_res9,antigen_res100)
        # interfaceRES9
        # interfaceRES100
    )





if __name__ == "__main__":

    start_time_readable = datetime.now()
    df_sabdab_all = get_sabdab_details()


    ab_both_chains = set(df_sabdab_all.query("Hchain == Hchain and Lchain == Lchain").pdb.values)
    print(len(ab_both_chains))
    # prots_fullab = prots.query(f"idcode.isin({list(ab_both_chains)})")

    
    base = "../structures/raw_new"
    f_epitopes = "../data/cdr_epitope.pickle"
    df = pd.read_pickle(f_epitopes)
    new_df = []
    c = 0
    c_noFail = 0
    structures = 0 #AB AG complexes
#############for the area test..

    # doneAntigens = set() 
    # doneAntibody = set() 
    # # areas=[]
    # nAntigens = 1 
    # nAntibody = 1
    currentABAG_complex = set()

#############
    # structures_noFail = 0

    currentStructure = ''
    currentStrNumb = 0
    start_time= time()
    end_time =0
    times = []
    # pdb2pqrFailure = False
    fp = open("failureList.txt", 'w')
    currentNumberProcessedPQRS =0
    totalNumberPQRs = 0
    noContact=0
    nComplexes=0
    processedCompexes = set()
    nPQRfailures=0
    missingBothCounter=0
    
    for prot in iterate_prots(df): #runs on cdr 1-2-3 of left chain then cdr 1-2-3 right chain
        #first structure is 7mhy to use for debug
        # if prot["idcode"] == "6kn9":
        #     continue
        
        # if (prot["idcode"] != '1afv') and (prot["idcode"] != '1ahw') and (prot["idcode"] != '1bj1')and (prot["idcode"] != '5lbs' ):
        #     continue
        # if (prot["idcode"] != '7mhy'):
        #     continue
        # if (prot["idcode"] != '6vlr'):
        #     continue

#NEW-----------------
        if(prot["idcode"]in ab_both_chains):
            pass
        else:
            if(currentStructure!=prot["idcode"]):
                print("Missing both chains..SKIPPING ", prot["idcode"])
                missingBothCounter +=1
                currentStructure=prot["idcode"]
            continue
    #############################
    #    
        # print(c) #internal counter over subsystems
        
        #DEBUG
        # if(prot["idcode"] != "3ogo"):
        #     continue
        ########
        if(currentStructure!=prot["idcode"]):
            currentNumberProcessedPQRS=0
       

        ab_chain_type = prot["chain_type"]
        if prot["chain_type"] == "K":
            ab_chain_type = "L"

        df_sabdab_prot = df_sabdab_all.query(
            f"pdb == '{prot['idcode']}' and {ab_chain_type}chain == '{prot['chainID']}'"
        )
        # try:
        if(len(df_sabdab_prot)!=0):
            antigen_type = df_sabdab_prot["antigen_type"].values[0]
        else:
            print("NOT IN DATABASE (90% SAAB)")
            continue
        # except:
        #     print(df_sabdab_prot["antigen_type"])
        #     print(prot['idcode'],ab_chain_type,prot['chainID'])
        #     input("error")
        antigen_chains = df_sabdab_prot["antigen_chain"].values[0].replace(" | ", "")
        if "protein" not in antigen_type:
            print('IGNORING:', antigen_type)
            continue

        

        antibody_chains = df_sabdab_prot[["Hchain", "Lchain"]].values[0]

        if prot["idcode"] == "4cc8" and antigen_chains == "G":
            print("IGNORING 4cc8")
            continue

        prot_dir = f'{EXPOSED_DIR}/{prot["idcode"]}/'
        if not os.path.isdir(prot_dir):
            os.mkdir(prot_dir)
        


        PROT_DIR = f'{EXPOSED_DIR}{prot["idcode"]}/'
        if not os.path.isdir(PROT_DIR):
            os.mkdir(PROT_DIR)

        (
            ab_ag_complex,
            cdrchain_ag_complex,
            cdr_ag_complex,
            antibody,
            antigen,
            cdrchain,cdr
        ) = get_systems(
            prot["idcode"],
            prot["cdr_atoms"],
            antigen_chains,
            prot["chainID"],
            antibody_chains,
            prot["cdr"],
        )
        
        if antibody_chains[0] != antibody_chains[0] or antibody_chains[1] != antibody_chains[1]:
            chainName =antibody_chains[0] if antibody_chains[1] != antibody_chains[1] else antibody_chains[1]
        else:
            chainName = ''.join(antibody_chains)
        
        complexName = f"{prot['idcode']}_complex_{chainName}_{antigen_chains}"#f"{idcode}_complex_{chains}_{antigen_chains}"

        # if(ab_ag_complex.fname in processedCompexes ):
        #         pass
        #     else:
        #         print("PDB2PQR fail: "+ complexName)
        #         failure = "PDB2PQR fail: "+ complexName
        #         fp.write('\n'+failure)
        #         processpedCompexes.add(ab_ag_complex.fname)
        #         nPQRfailures+=1
        #      #C
        #     continue
       
        if(ab_ag_complex is None):
            # pdb2pqrFailure = True
            # nameChain = "".join(antibody_chains)
            c+=1# KEEP COUNTING ALL SUBSYSTEM CONSIDERED
            if(complexName in processedCompexes ):
                pass
            else:
                processedCompexes.add(complexName)
                print("PDB2PQR fail: "+ complexName)
                failure = "PDB2PQR fail: "+ complexName
                fp.write('\n'+failure)  
                nPQRfailures+=1
             #C
            continue
        
        else:
            if(ab_ag_complex.fname in processedCompexes ):
                pass
            else:
                currentNumberProcessedPQRS +=1
                totalNumberPQRs+=1
                processedCompexes.add(ab_ag_complex.fname)
                # print("pqr name:",ab_ag_complex.fname)
                # print("current number of processed PQRs in the pdb:",currentNumberProcessedPQRS)
        
        #Qui crea il main file xyzr
        atom_res = pqr2xyzr(ab_ag_complex.fpqr, ab_ag_complex.fxyzr, prot["cdr"])
        
        run_nanoshaper(
            [ab_ag_complex, cdrchain_ag_complex, cdr_ag_complex, antibody, antigen,cdrchain,cdr]
        )

        (   failure,
            (ag_ab_interface, ag_ab_interface_res,ag_ab_interface_AREA),
            (ag_cdrchain_interface, ag_cdrchain_interface_res,ag_cdrchain_AREA),
            (ag_cdr_interface, ag_cdr_interface_res,ag_cdr_AREA),
            (ab_ag_interface, ab_ag_interface_res),
            (ag_interface, ag_interface_res,ag_interface_res9,ag_interface_res100)
        ) = get_interfaces(atom_res)

        c+=1
        # doneAntigens.add(antigen.fname)
        # doneAntibody.add(antibody.fname)
        # if(nAntigens!=len(doneAntigens)):
        #     nAntigens=len(doneAntigens)
        #     # av_antigenArea+=antigen.area
        #     # areas.append(antigen.area)
        #     # print(antigen.fname)
        #     # print("Area antigen = ",antigen.area)
       

        prot["ag_ab_interface"] = ag_ab_interface
        prot["ag_cdrchain_interface"] = ag_cdrchain_interface
        prot["ag_cdr_interface"] = ag_cdr_interface
        prot["ab_ag_interface"] = ab_ag_interface

        prot["ag_ab_interface_res"] = ag_ab_interface_res
        prot["ag_cdrchain_interface_res"] = ag_cdrchain_interface_res
        prot["ag_cdr_interface_res"] = ag_cdr_interface_res
        prot["ab_ag_interface_res"] = ab_ag_interface_res

        #NEW 
        prot["ag_ab_interface_AREA"] = ag_ab_interface_AREA #0.5*(A_antigen - (A_complex-A_antibody))
        # prot["ag_ab_interface_AREA"] = ab_ag_interface_AREA
        prot["ag_cdrchain_AREA"] = ag_cdrchain_AREA
        prot["ag_cdr_AREA"] = ag_cdr_AREA

        # prot["antigenName"] = f"{antigen_chains}" #--> to orient oneself better
        prot["contact name"] = antibody.fname+'_'+antigen.chainName
        prot["ag_interface"] = ag_interface
        prot["ag_interface_res"] = ag_interface_res
        prot["ag_interface_res_rp9"] = ag_interface_res9
        prot["ag_interface_res_rp100"] = ag_interface_res100



        new_df.append(prot.to_dict())

        end_time = time()
        if(ab_ag_complex.fname in currentABAG_complex):
            pass
        else:
            currentABAG_complex.add(ab_ag_complex.fname)
            print(ab_ag_complex.fname)
            print("Interface AREA = ", ag_ab_interface_AREA)

        
        # input('continue?')
        # print(nAntigens,len(doneAntigens))
        # if(len(doneAntigens)==12):
        #     av_antigenArea = np.average(areas)
        #     print(areas)
        #     print(av_antigenArea)
        #     break


        if(structures==12):
            break
        

        if(failure is not None):
            if(currentStrNumb!=structures):
                currentStrNumb=structures
                fp.write("\n")
            fp.write(failure+" ")
            if(failure=='NO AB AG CONTACT!'+ab_ag_complex.fname):
                noContact+=1 
                # NOTE: does NOT include PQR failures
            else:
                pass
                # nComplexes+=1
        else:
            c_noFail += 1
            # nComplexes +=1

        if(currentStructure!=prot["idcode"]):
            print(prot['idcode'])
            structures+=1
            
            # print("Number processed complexes in the structure =",currentNumberProcessedPQRS)
            print("Total number of processed PQRs",totalNumberPQRs)
            print("Total number of complexes not contact",noContact)
            # print("Total number of complexes ",nComplexes)
            
            # currentNumberProcessedPQRS=0
            # structures_noFail+=1
            
            # print(ab_ag_complex.fname)
            print("structure Number: ",structures)
            currentStructure = prot["idcode"] #complexName#prot["idcode"]
            
            
            elapsed = end_time-start_time
            times.append(elapsed)
            print("ELAPSED TIME= ",elapsed)
            start_time = time()
            # print(c,c_noFail)
        
        
    # fp.write(repr(failureList).replace("[","\n").replace("]","").replace(", ","\n"))
    fp.close()
    print("TOTAL NUMBER OF SYSTEMS CONSIDERED (including all sub-systems): ",c)
    print("TOTAL NUMBER OF PDB CONSIDERED: ",structures)
    print("TOTAL NUMBER OF PROCESSED SYSTEMS (including all sub-systems) WITH CONTACT: ",c_noFail)
    print()
    print("NUMBER OF SUCCESFULY PROCESSED PQRS",totalNumberPQRs)
    print("NUMBER OF FAILED PQRS",nPQRfailures)
    print("NUMBER OF COMPLEXES NOT IN CONTACT ",noContact)
    print("NUMBER OF SKIPPED PDBS SINCE NO BOTH CHAINS ",missingBothCounter)
    
    # print("NUMBER OF PROCESSED IN CONTACT ",nComplexes)
    
    # print("TOTAL NUMBER OF SUCCESSFULLY PROCESSED AB-AG COMPLEXES: ",structures_noFail)
    print("TOTAL RUN TIME = %.2f seconds" %(sum(times)))
    end_time_readable=datetime.now()
    print('Duration: {}'.format(end_time_readable - start_time_readable))

    df = pd.DataFrame(new_df)
    df.to_pickle("../data/epitope_buried_withAREA.pickle")
