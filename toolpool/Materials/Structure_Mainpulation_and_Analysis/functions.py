import pickle
from typing import List, Dict, Union, Any
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def save_to_pickle(data: Any, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Save data to a pickle file.
    Parameters:
        data: Any, the data to be saved.
        file_path: str, the path to the output pickle file where the data will be saved.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> Any:
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

def create_structure_from_file(structure_file_path: str, output_structure_pickle_file_path: str) -> None:
    """
    Name: create_structure_from_file
    Description: Creates a structure object from a given file and saves it to a pickle file.
    Parameters:
        structure_file_path: str, path to the structure file (e.g., CIF, POSCAR).
        output_structure_pickle_file_path: str, path to the output pickle file where the structure object will be saved.
    Returns:
        None
    """
    structure = Structure.from_file(structure_file_path)
    save_to_pickle(structure, output_structure_pickle_file_path)

def create_structure_from_spacegroup(spacegroup: str, lattice: List[List[float]], species: List[str], coords: List[List[float]], output_structure_pickle_file_path: str) -> None:
    """
    Name: create_structure_from_spacegroup
    Description: Creates a structure object from spacegroup information and saves it to a pickle file.
    Parameters:
        spacegroup: str, the spacegroup symbol (e.g., "Fm-3m").
        lattice: List[List[float]], the lattice parameters (e.g., [[a, 0, 0], [0, b, 0], [0, 0, c]]).
        species: List[str], a list of species (e.g., ["Na", "Cl"]).
        coords: List[List[float]], a list of fractional coordinates (e.g., [[0, 0, 0], [0.5, 0.5, 0.5]]).
        output_structure_pickle_file_path: str, path to the output pickle file where the structure object will be saved.
    Returns:
        None
    """
    structure = Structure.from_spacegroup(spacegroup, lattice, species, coords)
    save_to_pickle(structure, output_structure_pickle_file_path)

def add_site_property(structure_pickle_file_path: str, property_name: str, values: List[Union[int, float, str]], output_structure_pickle_file_path: str) -> None:
    """
    Name: add_site_property
    Description: Adds a site property to the structure and saves the modified structure to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to the input pickle file containing the structure object.
        property_name: str, the name of the site property.
        values: List[Union[int, float, str]], a list of property values corresponding to each site.
        output_structure_pickle_file_path: str, path to the output pickle file where the modified structure will be saved.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    structure.add_site_property(property_name, values)
    save_to_pickle(structure, output_structure_pickle_file_path)

def remove_sites(structure_pickle_file_path: str, indices: List[int], output_structure_pickle_file_path: str) -> None:
    """
    Name: remove_sites
    Description: Removes sites from the structure based on given indices and saves the modified structure to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to the input pickle file containing the structure object.
        indices: List[int], a list of indices of sites to be removed.
        output_structure_pickle_file_path: str, path to the output pickle file where the modified structure will be saved.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    structure.remove_sites(indices)
    save_to_pickle(structure, output_structure_pickle_file_path)

def get_symmetry_dataset(structure_pickle_file_path: str, output_structure_pickle_file_path: str) -> None:
    """
    Name: get_symmetry_dataset
    Description: Retrieves the symmetry dataset for the given structure and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to the input pickle file containing the structure object.
        output_structure_pickle_file_path: str, path to the output pickle file where the symmetry dataset will be saved.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    analyzer = SpacegroupAnalyzer(structure)
    symmetry_dataset = analyzer.get_symmetry_dataset()
    save_to_pickle(symmetry_dataset, output_structure_pickle_file_path)

def get_primitive_structure(structure_pickle_file_path: str, output_structure_pickle_file_path: str) -> None:
    """
    Name: get_primitive_structure
    Description: Gets the primitive structure of the given structure and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to the input pickle file containing the structure object.
        output_structure_pickle_file_path: str, path to the output pickle file where the primitive structure will be saved.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    analyzer = SpacegroupAnalyzer(structure)
    primitive_structure = analyzer.get_primitive_structure()
    save_to_pickle(primitive_structure, output_structure_pickle_file_path)
