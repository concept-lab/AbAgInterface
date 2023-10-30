from pandas import concat
import pymeshlab as ml
ms = ml.MeshSet()

# ms.load_new_mesh('7mhy_antigenPatch.off')

# ms.generate_splitting_by_connected_components(delete_source_mesh=True)

# n_connected = ms.number_meshes() #number of connected meshes

# print(" Number of connected components = ",n_connected)
# # Area of each connected component:

# for id in range(1,n_connected+1):
#     print("meshID=",id)
#     ms.set_current_mesh(id)
#     area = ms.get_geometric_measures()['surface_area']
#     print('area = ', area)

thresholdFractionArea = 2#%

import os
import glob
import numpy as  np

print("%d %% fraction of patches will be discarded" %thresholdFractionArea)

dataDIR = '../../structures/reducedResults/'
# structures = [s[0] for s in os.walk(dataDIR)]
structures = [ f.name for f in os.scandir(dataDIR) if f.is_dir() ]
print(structures)
# largerAreas = []
# smallerAreas =[]
# connectedComponents = []
outDict={}
for pdb in structures:
    path = dataDIR+pdb+'/'
    ag_ab_contacts = [f for f in glob.glob(path+'*_AGcontactPatch.off')]
    #esempio 6xcj_antibody_HL_antigen_G_AGcontactPatch.off
    # matchLine = "(\d)"
    # ab_agName = re.match(matchLine,ag_ab_contacts)
    
    for contact in ag_ab_contacts:
        print('\n',contact)
        ab_agName= pdb+"_"+"_".join(contact.split('_')[1:5])
        # print(ab_agName)
        # input()
        ms.clear()
        ms.load_new_mesh(contact)
        ms.generate_splitting_by_connected_components(delete_source_mesh=True)
        n_connected = ms.number_meshes() #number of connected meshes
        # connectedComponents.append(n_connected)
        print(" Number of connected components = ",n_connected)
        areas=[]
        for id in range(1,n_connected+1):
            print("meshID=",id)
            ms.set_current_mesh(id)
            area = ms.get_geometric_measures()['surface_area']
            print('area = ', area)
            areas.append(area)
        areas = np.array(areas)
        total_area = np.sum(areas)
        fraction_areas = areas/total_area *100
        print("fraction areas:", fraction_areas)
        indx = (fraction_areas > thresholdFractionArea).nonzero()
        areas = areas[indx]
        # fraction_areas = fraction_areas[indx]
        #UPDATE IF NECESSAY
        if(len(indx)<len(fraction_areas)):
            print("\nDISCARDING TOO SMALL DISCONNECTED PATCH")
            total_area = np.sum(areas)
            fraction_areas = areas/total_area *100 
            n_connected = len(fraction_areas)
            print("--> Number of connected components = ",n_connected)
        # largerAreas.append(np.amax(areas))
        # smallerAreas.append(np.amin(areas))
        
        # outDict[ab_agName]={'n connected':n_connected,'areas':areas,'max_area':np.amax(areas),'min_area':np.amin(areas),'fraction total':fraction_areas}
        outDict[ab_agName]={'n connected':n_connected,'areas':areas,'fraction total':fraction_areas}
print(outDict)
        

# print(outDict)

#save pickle TODO
import pandas as pd
df = pd.DataFrame(outDict)
df.to_pickle("contacts.pickle")