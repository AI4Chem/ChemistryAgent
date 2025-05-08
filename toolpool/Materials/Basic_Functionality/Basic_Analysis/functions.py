import pickle
import pymatgen.core as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

def save_to_pickle(data: any, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Save data to a pickle file.
    Parameters:
        data: Any, the data to be saved.
        file_path: str, the path to the output pickle file where the data will be saved.
    Returns:
        None
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

def analyze_symmetry(structure_pickle_file_path: str) -> str:
    """
    Name: analyze_symmetry
    Description: Analyzes the symmetry of a given crystal structure and saves its space group symbol to a pickle file.
    Parameters:
        structure_pickle_file_path: str, the path to the input pickle file containing the structure object.
    Returns:
        str, the space group symbol of symmetry analyzed.
    """
    structure = load_from_pickle(structure_pickle_file_path)
    finder = SpacegroupAnalyzer(structure)
    space_group_symbol = finder.get_space_group_symbol()
    return space_group_symbol

def match_structures(structure_pickle_file_path_1: str, structure_pickle_file_path_2: str) -> bool:
    """
    Name: match_structures
    Description: Compares two structures to determine if they are topologically identical.
    Parameters:
        structure_pickle_file_path_1: str, path to the first input pickle file containing the structure object.
        structure_pickle_file_path_2: str, path to the second input pickle file containing the structure object.
    Returns:
        bool, whether two structures are topologically identical.
    """
    structure1 = load_from_pickle(structure_pickle_file_path_1)
    structure2 = load_from_pickle(structure_pickle_file_path_2)
    matcher = StructureMatcher()
    result = matcher.fit_anonymous(structure1, structure2)
    return result