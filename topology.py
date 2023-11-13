from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGSLineSearch
import numpy as np
import os
import sys
import warnings
warnings.filterwarnings("ignore")

def get_contact_matrix(atoms, cutoff):
    number_of_atoms = atoms.get_number_of_atoms()
    rij = atoms.get_all_distances(mic = True)
    contact_matrix = np.where(rij <= cutoff, 1, 0)
    edge_list = np.array(np.triu_indices(number_of_atoms, k = 1)).transpose()
    within_cutoff = np.where(contact_matrix[edge_list[:, 0], edge_list[:, 1]] == 1)[0]
    edge_list = edge_list[within_cutoff]
    return contact_matrix, edge_list

def get_molecule_set(edge_list, number_of_atoms, molecule_set = None, atom_belong_to_molecule = None):
    if atom_belong_to_molecule is None: atom_belong_to_molecule = [-1] * number_of_atoms
    if molecule_set is None: molecule_set = {}
    new_molecule_name = max(atom_belong_to_molecule) + 1
    for edge in edge_list:
        atom_0 = edge[0]
        atom_1 = edge[1]
        atom_0_molecule = atom_belong_to_molecule[atom_0]
        atom_1_molecule = atom_belong_to_molecule[atom_1]
        if atom_0_molecule == -1 and atom_1_molecule == -1:
            molecule_set[new_molecule_name] = [atom_0, atom_1]
            atom_belong_to_molecule[atom_0] = new_molecule_name
            atom_belong_to_molecule[atom_1] = new_molecule_name
            new_molecule_name += 1
        elif atom_0_molecule == -1:
            molecule_set[atom_1_molecule].append(atom_0)
            atom_belong_to_molecule[atom_0] = atom_1_molecule
        elif atom_1_molecule == -1:
            molecule_set[atom_0_molecule].append(atom_1)
            atom_belong_to_molecule[atom_1] = atom_0_molecule
        elif atom_0_molecule != atom_1_molecule:
            min_molecule_name = min(atom_0_molecule, atom_1_molecule)
            max_molecule_name = max(atom_0_molecule, atom_1_molecule)
            tmp_molecule = molecule_set.pop(max_molecule_name)
            for j in tmp_molecule: atom_belong_to_molecule[j] = min_molecule_name
            molecule_set[min_molecule_name] += tmp_molecule
    return molecule_set, atom_belong_to_molecule

cutoff = 1.2
traj_dir = int(sys.argv[1])
output_dir = int(sys.argv[2])
filename_list = os.popen(f"ls {traj_dir}/*.extxyz").read().split('\n')
os.system(f"mkdir -p {output_dir}")
filename_list.pop()
for filename in filename_list:

    output_filename = filename[filename.rindex("/") + 1:filename.rindex(".extxyz")] + "-cluster.extxyz"
    if os.path.isfile(output_filename) == True: continue
    
    print(filename, flush = True)

    data = read(filename, format = "extxyz", index = ":")
    number_of_atoms = data[0].get_number_of_atoms()

    # create ref contact matrix, edge list, molecule set
    ref_contact_matrix, ref_edge_list = get_contact_matrix(data[0], cutoff)
    total_molecule_set, total_atom_belong_to_molecule = get_molecule_set(ref_edge_list, number_of_atoms)
    ref_molecule_set = {v: total_molecule_set[v][:] for v in total_molecule_set.keys()}
    ref_atom_belong_to_molecule = total_atom_belong_to_molecule[:]
    previous_contact_matrix = np.copy(ref_contact_matrix)

    # create final contact matrix, edge list, and update ref molecule set from final molecule set
    contact_matrix, edge_list = get_contact_matrix(data[-1], cutoff)
    final_molecule_set = {v: ref_molecule_set[v][:] for v in ref_molecule_set.keys()}
    final_atom_belong_to_molecule = ref_atom_belong_to_molecule[:]
    final_molecule_set, final_atom_belong_to_molecule = get_molecule_set(edge_list, number_of_atoms, final_molecule_set, final_atom_belong_to_molecule)
    number_of_contact_difference = np.sum((contact_matrix - previous_contact_matrix) ** 2)/2
    
    record_key = []
    if number_of_contact_difference > 0:
        print(number_of_contact_difference)

        for key in ref_molecule_set.keys():
            if final_molecule_set.get(key) is not None and ref_molecule_set.get(key) != final_molecule_set.get(key):
                record_key.append(key)
                print(f"{key} {ref_molecule_set.get(key)} {final_molecule_set.get(key)}")
                atom_indices = final_molecule_set.get(key)
                chemical_symbols = data[0].get_chemical_symbols()
                # chemical_symbols = [chemical_symbols[v] for v in atom_indices]
                cell = data[0].get_cell()
                pbc = data[0].get_pbc()
                new_data = []
                for atoms in data:
                    new_atoms = Atoms(chemical_symbols, positions = atoms.get_positions(), cell = cell, pbc = pbc)
                    reaction = np.zeros(new_atoms.get_number_of_atoms())
                    reaction[atom_indices] = 1
                    new_atoms.set_array('cluster', reaction)
                    # new_positions = new_atoms.get_distances(0, np.arange(new_atoms.get_number_of_atoms()), mic = True, vector = True)
                    # new_atoms.set_positions(new_positions)
                    # new_atoms.translate((new_atoms.get_cell()/2.).diagonal())
                    # new_atoms.translate(-new_atoms.get_center_of_mass() + (new_atoms.get_cell()/2.).diagonal())
                    new_data.append(new_atoms)
                output_filename = filename[filename.rindex("/") + 1:filename.rindex(".extxyz")] + f"-cluster-{key}.extxyz"
                write(f"{output_dir}/{output_filename}", new_data, format = 'extxyz')

        changed = False
        for iframe, atoms in enumerate(data[1:]):
            contact_matrix, edge_list = get_contact_matrix(atoms, cutoff)
            number_of_contact_difference = np.sum((contact_matrix - previous_contact_matrix) ** 2)/2
            if number_of_contact_difference > 0:
                changed = True
                print(iframe, number_of_contact_difference)
                total_molecule_set, total_atom_belong_to_molecule = get_molecule_set(edge_list, number_of_atoms, total_molecule_set, total_atom_belong_to_molecule)
            previous_contact_matrix = np.copy(contact_matrix)
        if not changed:
            print("No topology change")
    else:
        print("No topology change")

    cluster_index = [0] * number_of_atoms
    i_cluster = 1
    for key in total_molecule_set.keys():
        if ref_molecule_set.get(key) != total_molecule_set.get(key):
            for j in total_molecule_set[key]: cluster_index[j] = i_cluster
            i_cluster += 1

        if key in record_key:
            print(f"{key} {ref_molecule_set.get(key)} {total_molecule_set.get(key)}")

    cluster_index = np.array(cluster_index).astype(np.int)
    
    for atoms in data:
        atoms.set_array('cluster', cluster_index)

    output_filename = filename[filename.rindex("/") + 1:filename.rindex(".extxyz")] + "-cluster.extxyz"
    write(f"{output_dir}/{output_filename}", data, format = 'extxyz')
