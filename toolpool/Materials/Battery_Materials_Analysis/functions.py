import pickle
from pymatgen.core.structure import Structure
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.apps.battery.plotter import VoltageProfilePlotter
from pymatgen.analysis.diffusion.neb.pathfinder import NEBPathfinder
# from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer
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

def screen_electrode_materials(structure_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: screen_electrode_materials
    Description: Screens a given structure for its potential as an electrode material.
    Parameters:
        structure_pickle_file_path: str, the path to a pickle file containing a structure object.
    Returns:
        dict, a dictionary containing the electrode screening results.
    """

    structure = load_from_pickle(structure_pickle_file_path)
    
    electrode = InsertionElectrode(structure)
    
    return {
        'electrode': electrode,
        'max_voltage': electrode.get_max_voltage(),
        'min_voltage': electrode.get_min_voltage(),
        'average_voltage': electrode.get_average_voltage()
    }

def predict_voltage_profile(entries_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: predict_voltage_profile
    Description: Predicts the voltage profile for a given set of entries.
    Parameters:
        entries_pickle_file_path: str, the path to a pickle file containing a list of computed entries.
    Returns:
        dict, a dictionary containing the voltage profile.
    """
    voltage_analyzer = VoltageProfilePlotter(entries_pickle_file_path)

    return {
        'voltage_profile': voltage_analyzer.get_voltage_profile()
    }

def analyze_ion_diffusion_pathways(structure_pickle_file_path: str, min_slab_size: float, min_vacuum_size: float) -> Dict[str, Any]:
    """
    Name: analyze_ion_diffusion_pathways
    Description: Analyzes ion diffusion pathways in a given structure.
    Parameters:
        structure_pickle_file_path: str, the structure path to the pickle file.
        min_slab_size: float, minimum size of the slab.
        min_vacuum_size: float, minimum size of the vacuum layer.
    Returns:
        dict, a dictionary containing the diffusion pathways and their properties.
    """

    structure = load_from_pickle(structure_pickle_file_path)
    
    pathfinder = NEBPathfinder(structure, min_slab_size, min_vacuum_size)
    diffusion_paths = pathfinder.get_pathways()
    
    return {
        'diffusion_paths': diffusion_paths,
        'number_of_paths': len(diffusion_paths)
    }
