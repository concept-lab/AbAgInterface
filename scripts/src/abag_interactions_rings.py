import numpy as np
import mdtraj as md
import itertools
import random
import logging
from collections import namedtuple
PiPiPair = namedtuple('PiPiPair', ['antibody', 'antigen'])
PionPair = namedtuple('PionPair', ['ring', 'ion'])
Ring = namedtuple('Ring', ['indices', 'serials', 'resSeq', 'resSeq_str',
                           'resname', 'chain_ID', 'chain_type', 'CDR'])
Atom = namedtuple('Atom', ['index', 'serial', 'element', 'is_sidechain', 'resSeq',
                  'resSeq_str', 'resname', 'chain_ID', 'chain_type', 'CDR'])

########################################
# Some useful data and parameters
ring_atoms = {
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']}

ring_atoms_pi_ion = {
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']}

cutoff_ring = .5
cutoff_angle_pipi = .85
cutoff_angle_pion = .7
########################################


def get_ions(interface_atoms_pdb):
    def ions(atoms_iterable, resnames, element):
        control_set = set()
        for atom in atoms_iterable:
            if (atom.element == element) and (atom.resname in resnames) and atom.is_sidechain:
                # Get only 1 ion per residue (ASP/GLU or LYS/ARG)
                unique_id = str(atom.resSeq) + atom.resname + atom.chain_ID
                if unique_id not in control_set:
                    control_set.add(unique_id)
                    yield atom

    ab_anions = []
    for ion in ions(interface_atoms_pdb.antibody.values(), {'ASP', 'GLU'}, 'O'):
        ab_anions.append(ion)
    ag_anions = []
    for ion in ions(interface_atoms_pdb.antigen.values(), {'ASP', 'GLU'}, 'O'):
        ag_anions.append(ion)
    ab_cations = []
    for ion in ions(interface_atoms_pdb.antibody.values(), {'LYS', 'ARG'}, 'N'):
        ab_cations.append(ion)
    ag_cations = []
    for ion in ions(interface_atoms_pdb.antigen.values(), {'LYS', 'ARG'}, 'N'):
        ag_cations.append(ion)

    # ab_anions = tuple([atom for atom in interface_atoms_pdb.antibody.values() if (
    #     atom.element == 'O') and (atom.resname in {'ASP', 'GLU'}) and atom.is_sidechain])

    # ag_anions = tuple([atom for atom in interface_atoms_pdb.antigen.values() if (
    #     atom.element == 'O') and (atom.resname in {'ASP', 'GLU'}) and atom.is_sidechain])

    # ab_cations = tuple([atom for atom in interface_atoms_pdb.antibody.values() if (
    #     atom.element == 'N') and (atom.resname in {'LYS', 'ARG'}) and atom.is_sidechain])

    # ag_cations = tuple([atom for atom in interface_atoms_pdb.antigen.values() if (
    #     atom.element == 'N') and (atom.resname in {'LYS', 'ARG'}) and atom.is_sidechain])

    return tuple(ab_anions), tuple(ag_anions), tuple(ab_cations), tuple(ag_cations)


def is_complete(r, ring_atoms):
    min_n_atoms = len(ring_atoms[r.name]) + 4
    return not (r.n_atoms < min_n_atoms)


def get_atomic_data_from_ring_residue(positions, residue, ring_atoms):
    CD2_xyz = np.empty(3)
    CG_xyz = np.empty(3)
    CoM = np.zeros(3)
    natoms_ring = 0
    for atom in residue.atoms:
        if atom.name in ring_atoms[residue.name]:
            CoM += positions[atom.index, :]
            natoms_ring += 1
            if atom.name == 'CG':
                CG_index = atom.index
                CG_xyz = positions[atom.index, :]
            elif atom.name == 'CD2':
                CD2_xyz = positions[atom.index, :]
    CoM_xyz = CoM / natoms_ring

    v_1 = CG_xyz - CoM_xyz
    v_2 = CD2_xyz - CoM_xyz
    normal_vtor_not_norm = np.cross(v_1, v_2)

    return CG_index, CoM_xyz, normal_vtor_not_norm / np.linalg.norm(
        normal_vtor_not_norm)


def get_pipi_residues_from_chain(positions, chain, ring_atoms):
    CG_rings = []
    CoM_pipis_xzy = []
    normal_vectors = []
    for residue in chain.residues:
        if residue.name in ring_atoms:
            if is_complete(residue, ring_atoms):
                CG_index, CoM_xyz, normal_vtor = get_atomic_data_from_ring_residue(
                    positions, residue, ring_atoms)
                CG_rings.append(CG_index)
                CoM_pipis_xzy.append(CoM_xyz)
                normal_vectors.append(normal_vtor)

    return CG_rings, CoM_pipis_xzy, normal_vectors


def get_ring_data(trj_in, ab_chains, ring_atoms):
    CG_rings_ab = []
    CG_rings_ag = []
    CoM_pipis_xyz_ab = []
    CoM_pipis_xyz_ag = []
    normal_vectors_ab = []
    normal_vectors_ag = []
    for chain in trj_in.topology.chains:
        if chain.chain_id in ab_chains:
            CG_rings_ab_chain, CoM_pipis_xyz_ab_chain, normal_vectors_ab_chain =\
                get_pipi_residues_from_chain(trj_in.xyz[0], chain, ring_atoms)

            CG_rings_ab.extend(CG_rings_ab_chain)
            CoM_pipis_xyz_ab.extend(CoM_pipis_xyz_ab_chain)
            normal_vectors_ab.extend(normal_vectors_ab_chain)
        else:
            CG_rings_ag_chain, CoM_pipis_xyz_ag_chain, normal_vectors_ag_chain =\
                get_pipi_residues_from_chain(trj_in.xyz[0], chain, ring_atoms)

            CG_rings_ag.extend(CG_rings_ag_chain)
            CoM_pipis_xyz_ag.extend(CoM_pipis_xyz_ag_chain)
            normal_vectors_ag.extend(normal_vectors_ag_chain)

    return ({"antibody": CG_rings_ab, "antigen": CG_rings_ag},
            {"antibody": CoM_pipis_xyz_ab, "antigen": CoM_pipis_xyz_ag},
            {"antibody": normal_vectors_ab, "antigen": normal_vectors_ag})


def get_pipi_interactions(
        trj_in, CG_rings, CoM_rings_xyz, normal_vectors, cutoff_distance=.5,
        cutoff_angle=.9):
    pipi_ring_pairs = []
    for i, (com_i, v_i) in enumerate(
        zip(CoM_rings_xyz["antibody"],
            normal_vectors["antibody"])):
        for j, (com_j, v_j) in enumerate(
            zip(CoM_rings_xyz["antigen"],
                normal_vectors["antigen"])):
            distance = np.linalg.norm(com_i - com_j)
            angle = np.abs(np.dot(v_i, v_j))
            if distance < cutoff_distance and angle > cutoff_angle:
                ring_i = trj_in.topology.atom(CG_rings["antibody"][i]).residue
                ring_j = trj_in.topology.atom(CG_rings["antigen"][j]).residue

                pipi_ring_pairs.append(
                    PiPiPair(antibody=ring_i, antigen=ring_j))
    return pipi_ring_pairs


def get_ion_ring_interactions(
    trj_in, CG_rings, ions, CoM_rings_xyz, normal_vectors,
        positions, cutoff_distance, cutoff_angle=.7):

    if not len(CG_rings) or not len(ions):
        logging.warning("No good (whole) THR/TYR/PHE residues or, more likely, "
                        "no anions or cations in one of the molecules.")
        return []
    ids_ions = [ion.index for ion in ions]

    ring_ion_pairs = np.array(
        list(itertools.product(CG_rings, ids_ions)))

    ring_ion_distancias = md.compute_distances(
        trj_in, ring_ion_pairs).reshape((len(CG_rings), len(ids_ions)))

    indices_close_ion_CG = np.where(ring_ion_distancias < cutoff_distance)
    ion_ring_pairs = []
    for i, j in zip(*indices_close_ion_CG):
        com_xyz = CoM_rings_xyz[i]
        normal = normal_vectors[i]
        ON_xyz = positions[ids_ions[j]]
        ON_vector = ON_xyz - com_xyz
        norm_ON_vector = ON_vector / np.linalg.norm(ON_vector)

        distance = np.linalg.norm(ON_xyz - com_xyz)
        angle = np.abs(np.dot(norm_ON_vector, normal))

        if distance < cutoff_distance and angle > cutoff_angle:
            ring = trj_in.topology.atom(CG_rings[i]).residue
            ion_ring_pairs.append(PionPair(ring=ring, ion=ions[j]))
    return ion_ring_pairs


def get_chain_info_from_residue(pdb_idcode, epitope_buried_cleaned, residue):
    # Framework residue:
    cdr = 0
    rows_for_this_chain = epitope_buried_cleaned.query(
        f"idcode == '{pdb_idcode}' and chainID == '{residue.chain.chain_id}'")
    for row in rows_for_this_chain.iterrows():
        try:
            cdr_begin = int(row[1].cdr_begin.strip())
        except ValueError:
            cdr_begin = int(row[1].cdr_begin.strip()[:-1])
        try:
            cdr_end = int(row[1].cdr_end.strip())
        except ValueError:
            cdr_end = int(row[1].cdr_end.strip()[:-1])

        if cdr_begin <= residue.resSeq <= cdr_end:
            cdr = row[1].cdr
            break

    chain_type = 'L' if rows_for_this_chain.chain_type.values[0] in (
        'K', 'L') else rows_for_this_chain.chain_type.values[0]
    assert chain_type in ('K', 'L', 'H'), f"{pdb_idcode} returned invalid "\
        f"chain type for ring residue: {residue}"

    return chain_type, cdr


def get_data_from_ring_ring(pdb_idcode, epitope_buried_cleaned, mdtraj_ring_ring):
    rings = []
    for par in mdtraj_ring_ring:

        chain_type, cdr = get_chain_info_from_residue(
            pdb_idcode, epitope_buried_cleaned, par.antibody)

        rings.append(PiPiPair(
            antibody=Ring(indices=tuple([atom.index for atom in par.antibody.atoms]),
                          serials=tuple([atom.serial for atom in par.antibody.atoms]),
                          resSeq=par.antibody.resSeq, resSeq_str="",
                          resname=par.antibody.name,
                          chain_ID=par.antibody.chain.chain_id,
                          chain_type=chain_type, CDR=cdr),

            antigen=Ring(indices=tuple([atom.index for atom in par.antigen.atoms]),
                         serials=tuple([atom.serial for atom in par.antigen.atoms]),
                         resSeq=par.antigen.resSeq, resSeq_str="",
                         resname=par.antigen.name,
                         chain_ID=par.antigen.chain.chain_id,
                         chain_type='.', CDR=-1)))

    return tuple(rings)


def get_data_from_ring_ab_ion_ag(pdb_idcode, epitope_buried_cleaned, mdtraj_ring_ion):

    pares = []
    for par in mdtraj_ring_ion:

        chain_type, cdr = get_chain_info_from_residue(
            pdb_idcode, epitope_buried_cleaned, par.ring)

        pares.append(PionPair(
            ring=Ring(indices=tuple([atom.index for atom in par.ring.atoms]),
                      serials=tuple([atom.serial for atom in par.ring.atoms]),
                      resSeq=par.ring.resSeq, resSeq_str="",
                      resname=par.ring.name,
                      chain_ID=par.ring.chain.chain_id,
                      chain_type=chain_type, CDR=cdr),

            ion=par.ion))

    return tuple(pares)


def get_data_from_ring_ag_ion_ab(pdb_idcode, epitope_buried_cleaned, mdtraj_ring_ion):

    pares = []
    for par in mdtraj_ring_ion:

        pares.append(PionPair(
            ring=Ring(indices=tuple([atom.index for atom in par.ring.atoms]),
                      serials=tuple([atom.serial for atom in par.ring.atoms]),
                      resSeq=par.ring.resSeq, resSeq_str="",
                      resname=par.ring.name,
                      chain_ID=par.ring.chain.chain_id,
                      chain_type='.', CDR=-1),

            ion=par.ion))

    return tuple(pares)


def draw_pi_rings(trj_in, df_dataset, pdb_idcode, ring_ion_pairs, filename):
    with open(f"tempo/{filename}", "w") as fil:
        fil.write(f'from pymol import cmd\n\n')
        fil.write(f'cmd.load("{pdb_idcode}.pdb")\n')

        for n, ring_ion in enumerate(ring_ion_pairs):
            linea = f''
            ring_ion_id = 'ring_ion_' + str(n + 1)
            fil.write(f'cmd.select("id ')
            for ri in ring_ion:
                ri_serial = trj_in.topology.atom(ri).serial
                linea += f'{ri_serial}+'
            fil.write(linea[0:-1])
            fil.write(f'")\n')
            fil.write(f'cmd.set_name("sele", "{ring_ion_id}")\n')
            fil.write(f'cmd.show("spheres", "{ring_ion_id}")\n')
            color = "%06x" % random.randint(0, 0xFFFFFF)
            fil.write(f'cmd.color("0x{color}", "{ring_ion_id}")\n')
