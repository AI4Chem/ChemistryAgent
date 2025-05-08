import pickle
from pymatgen.core import Structure,PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import PeriodicSite
from pymatgen.analysis.defects.core import Vacancy, Interstitial
from typing import List, Dict, Any, Union


def save_to_pickle(obj: Any, file_path: str) -> None:
    """
    Saves an object to a pickle file.
    
    Parameters:
        obj: Any - The object to be saved.
        file_path: str - The path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(obj, f)


def load_from_pickle(file_path: str) -> Any:
    """
    Loads an object from a pickle file.
    
    Parameters:
        file_path: str - The path to the pickle file.
    
    Returns:
        Any - The loaded object.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)


def create_structure(lattice: List[List[float]], species: List[str], coords: List[List[float]], output_structure_pickle_file_path: str) -> Structure:
    """
    Name: create_structure
    Description: Creates a crystal structure from lattice vectors, species, and coordinates. Then saves it to a pickle file.
    Parameters:
        lattice: List[List[float]], a 3x3 matrix representing the lattice vectors.
        species: List[str], a list of element species (e.g., ['Si', 'O']).
        coords: List[List[float]], a list of fractional coordinates corresponding to the species.
        output_structure_pickle_file_path: str, path to the pickle file to save the created structure.
    Returns:
        Structure, the created crystal structure.
    """
    structure = Structure(lattice, species, coords)
    save_to_pickle(structure, output_structure_pickle_file_path)
    return structure


def analyze_symmetry(structure_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: analyze_symmetry
    Description: Analyzes the symmetry of the given crystal structure.
    Parameters:
        structure_pickle_file_path: str, path to a pickle file containing the crystal structure.
    Returns:
        dict, a dictionary containing symmetry information such as space group and point group.
    """
    if isinstance(structure_pickle_file_path, str):
        structure = load_from_pickle(structure_pickle_file_path)
    
    analyzer = SpacegroupAnalyzer(structure)
    return {
        "space_group": analyzer.get_space_group_symbol(),
        "point_group": analyzer.get_point_group_symbol(),
        "lattice": structure.lattice
    }


def modify_site_occupancy(structure_pickle_file_path: str, site_index: int, occupancy: float, output_structure_pickle_file_path: str) -> Structure:
    """
    Name: modify_site_occupancy
    Description: Modifies the occupancy of a specified site in the crystal structure then saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to a pickle file containing the crystal structure.
        site_index: int, the index of the site to modify.
        occupancy: float, the new occupancy value (0 to 1).
        output_structure_pickle_file_path: str, path to the pickle file to save the modified structure.
    Returns:
        Structure, the modified crystal structure.
    """
    if isinstance(structure_pickle_file_path, str):
        structure = load_from_pickle(structure_pickle_file_path)

    if 0 <= occupancy <= 1:
        site = structure[site_index]
        modified_site = PeriodicSite(site.species, site.coords, site.lattice)
        modified_site.occupancy = occupancy
        structure[site_index] = modified_site
        save_to_pickle(structure, output_structure_pickle_file_path)
        return structure
    else:
        raise ValueError("Occupancy must be between 0 and 1.")


def create_vacancy(structure_pickle_file_path: str, site_index: int, output_vacancy_pickle_file_path: str) -> Vacancy:
    """
    Name: create_vacancy
    Description: Creates a vacancy at a specified site in the crystal structure and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to a pickle file containing the crystal structure where the vacancy will be created.
        site_index: int, the index of the site to create the vacancy.
        output_vacancy_pickle_file_path: str, path to the pickle file to save the created vacancy.
    Returns:
        Vacancy, the created vacancy instance.
    """
    if isinstance(structure_pickle_file_path, str):
        structure = load_from_pickle(structure_pickle_file_path)

    site = structure[site_index]
    vacancy = Vacancy(structure, site)
    save_to_pickle(vacancy, output_vacancy_pickle_file_path)
    return vacancy


def create_interstitial(structure_pickle_file_path: str, element_str: str, coords: List[float], output_interstitial_pickle_file_path: str) -> Interstitial:
    """
    Name: create_interstitial
    Description: Creates an interstitial defect in the crystal structure and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, the crystal structure where the interstitial will be created or path to the pickle file.
        element_str: str, the element of the interstitial.
        coords: List[float], the coordinates of the interstitial.
        output_interstitial_pickle_file_path: str, path to the pickle file to save the created interstitial.
    Returns:
        Interstitial, the created interstitial instance.
    """
    if isinstance(structure_pickle_file_path, str):
         structure = load_from_pickle(structure_pickle_file_path)
     
    interstitial = Interstitial(structure, element_str, coords)
    save_to_pickle(interstitial, output_interstitial_pickle_file_path)
    return interstitial
