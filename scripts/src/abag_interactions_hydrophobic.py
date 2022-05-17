import logging
import numpy as np
import itertools
import mdtraj as md
import networkx as nx
import string
import random
from collections import namedtuple
Atom = namedtuple('Atom', ['index', 'serial', 'element', 'is_sidechain', 'resSeq',
                  'resSeq_str', 'resname', 'chain_ID', 'chain_type', 'CDR'])

########################################
# Some useful data and parameters
cutoff_clusters = .45
cutoff_carbons = .45
########################################


def get_waters(topologia):

    ids_wat = []
    for at in topologia.atoms:
        if at.residue.name == 'HOH':
            ids_wat.append(at.index)

    return ids_wat


def is_shielded(
        positions, C_cdr_id, C_epi_id, surrounding_ON_ids, angle_cdr_ON=.85,
        angle_epi_ON=-.2):
    C_cdr_xyz = positions[C_cdr_id, :]
    C_epi_xyz = positions[C_epi_id, :]

    vec_C_C = C_epi_xyz - C_cdr_xyz
    n_vec_C_C = vec_C_C / np.linalg.norm(vec_C_C)

    for ON_id in surrounding_ON_ids:
        ON_xyz = positions[ON_id, :]
        vec_cdr_ON = ON_xyz - C_cdr_xyz
        n_vec_cdr_ON = vec_cdr_ON / np.linalg.norm(vec_cdr_ON)
        vec_epi_ON = ON_xyz - C_epi_xyz
        n_vec_epi_ON = vec_epi_ON / np.linalg.norm(vec_epi_ON)
        # Useful for debugging:
        dot_Cab_Cag_vs_Cab_ON = np.dot(n_vec_C_C, n_vec_cdr_ON)
        dot_Cab_Cag_vs_Cag_ON = np.dot(n_vec_C_C, n_vec_epi_ON)
        ###
        if ON_id == 257 and C_cdr_id == 4270:
            print(f" {C_epi_id=} -- {dot_Cab_Cag_vs_Cab_ON=} - {dot_Cab_Cag_vs_Cag_ON=}")
        ###
        if (dot_Cab_Cag_vs_Cab_ON > angle_cdr_ON) and (dot_Cab_Cag_vs_Cag_ON < angle_epi_ON):
            return True, ON_id

    return False, 0


def get_carbons_graph(trj_in, interface_atoms_pdb, cutoff):

    ab_carbons = [atom.index for atom in interface_atoms_pdb.antibody.values()
                  if atom.element == 'C']

    ag_carbons = [atom.index for atom in interface_atoms_pdb.antigen.values()
                  if atom.element == 'C']

    ab_polars = [atom.index for atom in interface_atoms_pdb.antibody.values()
                 if atom.element in ('N', 'O')]
    ag_polars = [atom.index for atom in interface_atoms_pdb.antigen.values()
                 if atom.element in ('N', 'O')]
    polars = ab_polars + ag_polars

    C_ON_pairs = np.array(list(itertools.product(ab_carbons, polars)))
    C_C_pairs = np.array(list(itertools.product(ab_carbons, ag_carbons)))

    C_ON_distancias = md.compute_distances(
        trj_in, C_ON_pairs).reshape((len(ab_carbons), len(polars)))

    C_C_distancias = md.compute_distances(
        trj_in, C_C_pairs).reshape((len(ab_carbons), len(ag_carbons)))

    G = nx.Graph()
    indices_close_C_C_distancias = np.where(C_C_distancias < cutoff)
    mask_close_C_ON_distancias = C_ON_distancias < cutoff
    shielding_atoms_serial = {}

    for i, j in zip(*indices_close_C_C_distancias):
        C_cdr_id = ab_carbons[i]
        C_epi_id = ag_carbons[j]
        surrounding_ON_ids = \
            [polars[i] for i in np.where(mask_close_C_ON_distancias[i, :])[0]]

        shielded, ON_id = is_shielded(trj_in.xyz[0], C_cdr_id, C_epi_id, surrounding_ON_ids)

        if shielded:
            shielding_atoms_serial[i] = trj_in.topology.atom(ON_id).serial
        else:
            Cab = trj_in.topology.atom(C_cdr_id).serial
            Cag = trj_in.topology.atom(C_epi_id).serial

            if ON_id == 1091 and C_cdr_id == 1470:
                print(f" {C_epi_id=} - {C_cdr_id=} - {ON_id=}")

            G.add_edge(interface_atoms_pdb.antibody[Cab], interface_atoms_pdb.antigen[Cag])

    return G


def get_putative_clusters(G):
    pre_clusteres = []
    for cluster in sorted(nx.connected_components(G), key=len, reverse=True):
        pre_clusteres.append(cluster)

    return pre_clusteres


def merge_clusters(trj_in, pre_clusteres, cutoff_clusters):
    def clusters_are_close(trj_in, cluster_1, cluster_2, cutoff_clusters):
        cluster_2_indices = [carbon.index for carbon in cluster_2]

        for carbon in cluster_1:
            distancias = md.compute_distances(trj_in, list(
                itertools.product([carbon.index], cluster_2_indices)))

            close_clusters = np.any(distancias < cutoff_clusters)
            if close_clusters:
                return True

        return False

    H = nx.Graph()
    H.add_node(0)
    for i, clu_i in enumerate(pre_clusteres):
        for j in range(i + 1, len(pre_clusteres)):
            if clusters_are_close(
                    trj_in, clu_i, pre_clusteres[j],
                    cutoff_clusters):
                H.add_edge(i, j)
            else:
                H.add_node(i)
                H.add_node(j)

    clusteres = []
    for connected_clusters in sorted(
            nx.connected_components(H),
            key=len, reverse=True):
        new_cluster = []
        for i in connected_clusters:
            new_cluster.extend(pre_clusteres[i])
        clusteres.append(new_cluster)

    # Make sure that the clusters are sorted by size
    idx = np.flip(np.argsort([len(c) for c in clusteres]))
    sorted_clusters = []
    for i in idx:
        sorted_clusters.append(tuple(clusteres[i]))

    return tuple(sorted_clusters)


def draw_clusters(pdb_filename, interface_atoms_pdb, ag_chains, clusters,
                  filename):
    with open(filename, "w") as fil:
        fil.write(f'from pymol import cmd\n\n')
        fil.write(f'cmd.set("sphere_scale", "0.9")\n\n')
        fil.write(f'cmd.load("{pdb_filename}")\n')
        fil.write(f'cmd.color("salmon", "')
        for chainID in ag_chains:
            if chainID != '.':
                fil.write(f'chain {chainID} or ')
        fil.write(f'chain {ag_chains[-1]}")\n')
        fil.write(f'cmd.color("atomic", "(not elem C)")\n\n')

        # Show epitope residues as lines
        epitope_residues = set([atom.resSeq_str
                                for atom in interface_atoms_pdb.antigen.values()])

        for resi in epitope_residues:
            fil.write(f'cmd.show("lines", "resi {resi} and chain {ag_chains[0]}")' + '\n')

        # Show paratope residues as lines
        paratope_residues = set([(atom.resSeq_str, atom.chain_ID)
                                 for atom in interface_atoms_pdb.antibody.values()])

        for resi, chainID in paratope_residues:
            fil.write(
                f'cmd.show("lines", "resi {resi} and chain {chainID}")\n')

        # Finally, show interacting carbons as spheres
        for n, cluster in enumerate(clusters):
            linea = f''
            cluster_id = 'cluster_' + str(n + 1)
            fil.write(f'cmd.select("id ')
            for c in cluster:
                linea += f'{c.serial}+'
            fil.write(linea[0:-1])
            fil.write(f'")\n')
            fil.write(f'cmd.set_name("sele", "{cluster_id}")\n')
            fil.write(f'cmd.show("spheres", "{cluster_id}")\n')
            color = "%06x" % random.randint(0, 0xFFFFFF)
            fil.write(f'cmd.color("0x{color}", "{cluster_id}")\n')
            # fil.write(f'cmd.color("Gray", "{cluster_id}")\n')
