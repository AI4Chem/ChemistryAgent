import pickle
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer
from pymatgen.entries.computed_entries import ComputedEntry
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

def get_phase_diagram_properties(entries_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: get_phase_diagram_properties
    Description: Generates a phase diagram and returns key properties like kinds of entries(stable/unstable) and hull.
    Parameters:
        entries_pickle_file_path: str, the path to the pickle file of computed entries.
    Returns:
        dict, a dictionary containing key properties of the phase diagram.
    """
    entries = load_from_pickle(entries_pickle_file_path)
    
    phase_diagram = PhaseDiagram(entries)
    
    return {
        'stable_entries': phase_diagram.stable_entries,
        'unstable_entries': phase_diagram.unstable_entries,
        'hull': phase_diagram.get_hull_energy_dict()
    }

def analyze_diffusion(entries_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: analyze_diffusion
    Description: Analyzes diffusion properties and returns key properties like diffusion coefficient and activation energy.
    Parameters:
        entries_pickle_file_path: str, the path to the pickle file of computed entries.
    Returns:
        dict, a dictionary containing key properties of the diffusion analysis.
    """
    if isinstance(entries_pickle_file_path, str):
        entries = load_from_pickle(entries_pickle_file_path)
    
    diffusion_analyzer = DiffusionAnalyzer(entries)

    return {
        'diffusion_coefficient': diffusion_analyzer.diffusivity,
        'activation_energy': diffusion_analyzer.activation_energy
    }

def calculate_magnetic_properties(entries_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: calculate_magnetic_properties
    Description: Calculates magnetic properties and returns key properties like total magnetization and magnetic moments.
    Parameters:
        entries_pickle_file_path: str, the path to the pickle file of computed entries.
    Returns:
        dict, a dictionary containing key properties of the magnetic analysis.
    """

    entries = load_from_pickle(entries_pickle_file_path)
    
    magnetic_properties = []
    for entry in entries:
        structure = entry.get("structure")  # Ensure this corresponds to your data format
        if isinstance(structure, Structure):
            analyzer = CollinearMagneticStructureAnalyzer(structure)
            magnetic_properties.append({
                "total_magnetization": analyzer.total_magnetization,
                "magnetic_moments": analyzer.magnetic_moments,
                "is_magnetic": analyzer.is_magnetic,
                "magnetic_ordering": analyzer.ordering.value
            })
    
    return {
        "magnetic_properties": magnetic_properties
    }
