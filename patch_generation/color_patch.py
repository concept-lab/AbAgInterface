import numpy as np
import open3d as o3d
import os
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import time

#7mhy_antigen_A
#7mhy_antibody_MN
input_folder = "../structures/exposed/7mhy/"
input_protein = "7mhy_antigen_A"
input_ligand = "7mhy_cdrcomplex_M3_A_vertMapResBased"

# palette_extended = np.random.randint(256, size=(50, 3))
palette_extended = np.empty((4,3))
palette_extended[0,:] = [128, 128, 128] #none
palette_extended[1,:] = [255,0,255] #all anitgen antibody interface VIOLET
# palette_extended[1,:] = [0,204,255] #all anitgen antibody interface BLUE
palette_extended[2,:] = [255,204,0] #CDR chain YELLOW
palette_extended[3,:] = [255,0,0] #CDR  RED

print("Reading the SES triangle mesh ...")
mesh = o3d.io.read_triangle_mesh(input_folder + input_protein + ".off")
V = np.asarray(mesh.vertices)
F = np.asarray(mesh.triangles)

# Vertices that potentially require to be coloured:
V_labels = []
with open(input_folder + input_ligand + ".txt") as f:
    for line in f:
        V_labels = V_labels + [int(i) for i in line.split()]
V_labels = np.array(V_labels)


output_file = "./triangulatedSurf"
with open(output_file + ".txt", 'w') as file:
    file.write("COFF\n")
    file.write(f"{V.shape[0]} {F.shape[0]} 0\n")
    for idx, V_row in enumerate(V):
        file.write(f"{V_row[0]} {V_row[1]} {V_row[2]} ")
        file.write(f"{int(palette_extended[V_labels[idx]][0])} ")
        file.write(f"{int(palette_extended[V_labels[idx]][1])} ")
        file.write(f"{int(palette_extended[V_labels[idx]][2])} ")
        file.write(f"255\n")
    for F_row in F:
        file.write(f"3 {F_row[0]} {F_row[1]} {F_row[2]}\n")
os.rename(output_file + ".txt", output_file + ".off")
