import pickle
import pymatgen.core as mg
from typing import List, Union

def save_to_pickle(data: any, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Save data to a pickle file.
    Parameters:
        data: Any, the data to be saved.
        file_path: str, the path to the output pickle file where the data will be saved.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> any:
    """
    Name: load_from_pickle
    Description: Load data from a pickle file.
    Parameters:
        file_path: str, the path to the input pickle file from which the data will be loaded.
    Returns:
        data: Any, the data loaded from the pickle file.
    """
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def create_cubic_lattice(lattice_parameter: float, output_lattice_pickle_file_path: str) -> None:
    """
    Name: create_cubic_lattice
    Description: Creates a cubic lattice with the specified lattice parameter and saves it to a pickle file.
    Parameters:
        lattice_parameter: float, the lattice parameter for the cubic lattice.
        output_lattice_pickle_file_path: str, the path to the output pickle file where the lattice will be saved.
    Returns:
        None
    """
    lattice = mg.Lattice.cubic(lattice_parameter)
    save_to_pickle(lattice, output_lattice_pickle_file_path)

def create_structure(lattice_pickle_file_path: str, species: List[str], coordinates: List[List[float]], output_structure_pickle_file_path: str) -> None:
    """
    Name: create_structure
    Description: Creates a structure from a given lattice, species, and coordinates, and saves it to a pickle file.
    Parameters:
        lattice_pickle_file_path: str, The path to the input pickle file containing the lattice.
        species: List[str], a list of species (elements) present in the structure.
        coordinates: List[List[float]], a list of fractional coordinates corresponding to the species.
        output_structure_pickle_file_path: str, the path to the output pickle file where the structure will be saved.
    Returns:
        None
    """
    lattice = load_from_pickle(lattice_pickle_file_path)
    structure = mg.Structure(lattice, species, coordinates)
    save_to_pickle(structure, output_structure_pickle_file_path)

def modify_structure(structure_pickle_file_path: str, output_structure_pickle_file_path: str, operation: str, target: list) -> None:
    """
    Name: modify_structure
    Description: Modifies a given structure based on the specified operation and saves the modified structure to a pickle file.
    Parameters:
        structure_pickle_file_path: str, The path to the input pickle file containing the structure.
        output_structure_pickle_file_path: str, The path to the output pickle file where the modified structure will be saved.
        operation: str, The modification operation (e.g., 'make_supercell', 'delete', 'append', 'change', 'shift').
        target: list, target required for the modification operation.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    
    if operation == 'make_supercell':
        structure.make_supercell(target[0])  # args[0] is expected to be a list of integers for scaling
    elif operation == 'delete':
        del structure[target[0]]  # args[0] is the index to delete
    elif operation == 'append':
        structure.append(target[0], target[1])  # args[0] is the species, args[1] is the coordinates
    elif operation == 'change':
        structure[-1] = target[0]  # args[0] is the new species
    elif operation == 'shift':
        structure[target[0]] = target[1]  # args[0] is the index, args[1] is the new coordinates
    
    save_to_pickle(structure, output_structure_pickle_file_path)

def create_immutable_structure(structure_pickle_file_path: str, output_structure_pickle_file_path: str) -> None:
    """
    Name: create_immutable_structure
    Description: Converts a mutable structure to an immutable structure and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, the path to the input pickle file containing the mutable structure.
        output_structure_pickle_file_path: str, the path to the output pickle file where the immutable structure will be saved.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    immutable_structure = mg.IStructure.from_sites(structure)
    save_to_pickle(immutable_structure, output_structure_pickle_file_path)