import pickle
import pymatgen.core as mg
from pymatgen.io.vasp.sets import MPRelaxSet

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

def write_structure_to_file(structure_pickle_file_path: str, output_structure_file_path: str) -> None:
    """
    Name: write_structure_to_file
    Description: Writes a structure to a specified file format. The format is intelligently determined from the file name extension if not provided.
    Parameters:
        structure_pickle_file_path: str, the path to the input pickle file containing the structure.
        output_structure_file_path: str, the name of the file to write the structure to.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    structure.to(filename=output_structure_file_path)

def read_structure_from_file(structure_file_path: str, structure_pickle_file_path: str) -> None:
    """
    Name: read_structure_from_file
    Description: Reads a structure from a specified file and saves it to a pickle file.
    Parameters:
        structure_file_path: str, the name of the file to read the structure from.
        structure_pickle_file_path: str, the path to the output pickle file where the structure will be saved.
    Returns:
        None
    """
    structure = mg.Structure.from_file(structure_file_path)
    save_to_pickle(structure, structure_pickle_file_path)

def create_vasp_input_files(structure_pickle_file_path: str, output_dir: str) -> None:
    """
    Name: create_vasp_input_files
    Description: Generates a complete set of VASP input files for the given structure.
    Parameters:
        structure_pickle_file_path: str, the path to the input pickle file containing the structure.
        output_dir: str, The directory where the input files will be written.
    Returns:
        None
    """
    structure = load_from_pickle(structure_pickle_file_path)
    vasp_input_set = MPRelaxSet(structure)
    vasp_input_set.write_input(output_dir)