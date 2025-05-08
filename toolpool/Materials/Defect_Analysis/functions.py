import pickle
from pymatgen.analysis.defects.core import Defect
from pymatgen.analysis.defects.thermo import DefectEntry
from pymatgen.core import Structure

def save_to_pickle(data, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Saves the data to a pickle file.
    Parameters:
        data: any - The data to be saved.
        file_path: str - The path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> any:
    """
    Name: load_from_pickle
    Description: Loads data from a pickle file.
    Parameters:
        file_path: str - The path to the pickle file.
    Returns:
        any - The loaded data.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def create_defect(structure_pickle_file_path: str, defect_type: str, defect_site: tuple, multiplicity: int, charge: int, output_defect_pickle_file_path: str) -> Defect:
    """
    Name: create_defect
    Description: Creates a defect instance from given parameters.
    Parameters:
        structure_pickle_file_path: str, path to the structure file.
        defect_type: str, type of the defect (e.g., 'vacancy', 'interstitial', 'substitution').
        defect_site: tuple, site of the defect in fractional coordinates.
        multiplicity: int, multiplicity of the defect.
        charge: int, charge state of the defect.
        output_defect_pickle_file_path: str, the path to a pickle file to save the Defect instance created.
    Returns:
        Defect, the created Defect instance.
    """
    structure = Structure.from_file(structure_pickle_file_path)
    defect = Defect(structure, defect_type, defect_site, multiplicity, charge)
    save_to_pickle(defect, output_defect_pickle_file_path)
    return defect

def load_defect_from_pickle(defect_pickle_file_path: str) -> Defect:
    """
    Name: load_defect_from_pickle
    Description: Loads a Defect instance from a pickle file.
    Parameters:
        defect_pickle_file_path: str, the path to the pickle file containing a Defect instance.
    Returns:
        Defect, the loaded Defect instance.
    """
    return load_from_pickle(defect_pickle_file_path)

def get_defect_type(defect_pickle_file_path: str) -> str:
    """
    Name: get_defect_type
    Description: Retrieves the type of the defect.
    Parameters:
        defect_pickle_file_path: str, the path to a pickle file containing a Defect instance.
    Returns:
        str, the defect type.
    """
    return load_defect_from_pickle(defect_pickle_file_path).defect_type

def get_defect_site(defect_pickle_file_path: str) -> tuple:
    """
    Name: get_defect_site
    Description: Retrieves the site of the defect in fractional coordinates.
    Parameters:
        defect_pickle_file_path: str, the path to a pickle file containing a Defect instance.
    Returns:
        tuple, the defect site.
    """
    return load_defect_from_pickle(defect_pickle_file_path).defect_site

def get_defect_charge(defect_pickle_file_path: str) -> int:
    """
    Name: get_defect_charge
    Description: Retrieves the charge state of the defect.
    Parameters:
        defect_pickle_file_path: str, the path to a pickle file containing a Defect instance.
    Returns:
        int, the defect charge state.
    """
    return load_defect_from_pickle(defect_pickle_file_path).charge

def create_defect_entry(defect_pickle_file_path: str, energy: float, output_defect_entry_pickle_file_path: str) -> DefectEntry:
    """
    Name: create_defect_entry
    Description: Creates a DefectEntry instance from given parameters and saves it to a pickle file.
    Parameters:
        defect_pickle_file_path: str, the path to the pickle file containing a Defect instance.
        energy: float, the energy associated with the defect.
        output_defect_entry_pickle_file_path: str, the path to a pickle file to save the created DefectEntry instance.
    Returns:
        DefectEntry, the created DefectEntry instance.
    """
    defectentry = DefectEntry(load_defect_from_pickle(defect_pickle_file_path), energy)
    save_to_pickle(defectentry, output_defect_entry_pickle_file_path)
    return defectentry

def load_defect_entry_from_pickle(defect_entry_pickle_file_path: str) -> DefectEntry:
    """
    Name: load_defect_entry_from_pickle
    Description: Loads a DefectEntry instance from a pickle file.
    Parameters:
        defect_entry_pickle_file_path: str, the path to a pickle file containing a DefectEntry instance.
    Returns:
        DefectEntry, the loaded DefectEntry instance.
    """
    return load_from_pickle(defect_entry_pickle_file_path)

def get_defect_energy(defect_entry_pickle_file_path: str) -> float:
    """
    Name: get_defect_energy
    Description: Retrieves the energy associated with the defect.
    Parameters:
        defect_entry_pickle_file_path: str, the path to a pickle file containing a DefectEntry instance.
    Returns:
        float, the defect energy.
    """
    return load_defect_entry_from_pickle(defect_entry_pickle_file_path).energy