import pickle
from pymatgen.core.composition import Composition
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


def analyze_elemental_composition(composition_str: str) -> Dict[str, float]:
    """
    Name: analyze_elemental_composition
    Description: Analyzes the elemental composition of a given formula string.
    Parameters:
        composition_str: str, the chemical formula as a string.
    Returns:
        dict, a dictionary with element symbols as keys and their amounts as values.
    """
    composition = Composition(composition_str)
    return composition.get_el_amt_dict()


def guess_oxidation_states(composition_str: str) -> List[Dict[str, int]]:
    """
    Name: guess_oxidation_states
    Description: Guesses the possible oxidation states for the given composition.
    Parameters:
        composition_str: str, the chemical formula as a string.
    Returns:
        list, a list of dictionaries with elements as keys and their guessed oxidation states as values.
    """
    composition = Composition(composition_str)
    return composition.oxi_state_guesses()


def mix_elements(composition_str: str, element_ratios: Dict[str, float], output_composition_pickle_file_path: str) -> Composition:
    """
    Name: mix_elements
    Description: Mixes elements in the given composition with specified ratios.
    Parameters:
        composition_str: str, the chemical formula as a string.
        element_ratios: dict, a dictionary with element symbols as keys and their fractional ratios as values.
        output_composition_pickle_file_path: str, path to the pickle file to save the created composition object.
    Returns:
        Composition, the mixed composition.
    """
    composition = Composition(composition_str)

    mixed_composition = composition.fractional_composition
    save_to_pickle(mixed_composition, output_composition_pickle_file_path)
    return mixed_composition
