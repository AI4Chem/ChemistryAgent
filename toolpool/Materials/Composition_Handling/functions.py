import pickle
from pymatgen.core import Composition
from typing import Dict

def create_composition(composition_str: str, output_composition_pickle_file_path: str) -> None:
    """
    Name: create_composition
    Description: Creates a Composition object from a given chemical formula and saves it to a pickle file.
    Parameters:
        composition_str: str, the chemical formula (e.g., "NaCl", "CaCO3").
        output_composition_pickle_file_path: str, the path to the output pickle file where the Composition object will be saved.
    Returns:
        None
    """
    composition = Composition(composition_str)
    with open(output_composition_pickle_file_path, 'wb') as f:
        pickle.dump(composition, f)

def get_atomic_fraction(composition_pickle_file_path: str) -> Dict:
    """
    Name: get_atomic_fraction
    Description: Retrieves the atomic fraction of each element in the composition.
    Parameters:
        composition_pickle_file_path: str, the path to the input pickle file containing the Composition object.
    Returns:
        Dict
    """
    with open(composition_pickle_file_path, 'rb') as f:
        composition = pickle.load(f)

    atomic_fractions: Dict[str, float] = {el.symbol: composition.get_atomic_fraction(el) for el in composition.elements}
    return atomic_fractions

def get_weight_fraction(composition_pickle_file_path: str) -> Dict:
    """
    Name: get_weight_fraction
    Description: Retrieves the weight fraction of each element in the composition.
    Parameters:
        composition_pickle_file_path: str, the path to the input pickle file containing the Composition object.
    Returns:
        Dict
    """
    with open(composition_pickle_file_path, 'rb') as f:
        composition = pickle.load(f)

    weight_fractions: Dict[str, float] = {el.symbol: composition.get_wt_fraction(el) for el in composition.elements}
    return weight_fractions

def get_reduced_composition(composition_pickle_file_path: str, output_composition_pickle_file_path: str) -> None:
    """
    Name: get_reduced_composition
    Description: Retrieves the reduced composition of the given composition and saves it to a pickle file.
    Parameters:
        composition_pickle_file_path: str, the path to the input pickle file containing the Composition object.
        output_composition_pickle_file_path: str, the path to the output pickle file where the reduced composition will be saved.
    Returns:
        None
    """
    with open(composition_pickle_file_path, 'rb') as f:
        composition = pickle.load(f)

    reduced_composition = composition.reduced_composition
    with open(output_composition_pickle_file_path, 'wb') as f:
        pickle.dump(reduced_composition, f)
