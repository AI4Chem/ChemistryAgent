import pickle
from pymatgen.io.cif import CifWriter, CifParser
from pymatgen.io.vasp import VaspInput, Vasprun
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure

def save_to_pickle(data, file_path: str) -> None:
    """
    Saves the data to a pickle file.
    
    Parameters:
        data: any - The data to be saved.
        file_path: str - The path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> any:
    """
    Loads data from a pickle file.
    
    Parameters:
        file_path: str - The path to the pickle file.
    
    Returns:
        any - The loaded data.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def handle_cif(structure_pickle_file_path: str, output_cif_pickle_file_path: str):
    """
    Name: handle_cif
    Description: Handles CIF file operations and saves it to a pickle file.
    Parameters:
        structure_pickle_file_path: str, path to the structure pickle file.
        output_cif_pickle_file_path: str, path to the CIF file to write to.
    Returns:
        Structure, the structure loaded or None if writing to CIF.
    """
    writer = CifWriter(load_from_pickle(structure_pickle_file_path))
    writer.write_file(output_cif_pickle_file_path)


def handle_vasp(vasp_pickle_file_path: str):
    """
    Name: handle_vasp
    Description: Handles VASP input operations.
    Parameters:
        vasp_pickle_file_path: str, path to the pickle file containing a saved VaspInput or Vasprun instance.
    Returns:
        VaspInput or Vasprun, the loaded instance.
    """
    return load_from_pickle(vasp_pickle_file_path)


def get_structure_by_material_id(api_key: str, material_id: str) -> Structure:
    """
    Name: get_structure_by_material_id
    Description: Retrieves a structure by material ID from the Materials Project.
    Parameters:
        api_key: str, the API key for accessing the Materials Project.
        material_id: str, the material ID.
    Returns:
        Structure, the retrieved structure.
    """
    with MPRester(api_key) as rester:
        return rester.get_structure_by_material_id(material_id)
