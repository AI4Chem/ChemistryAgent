import pickle
#from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
#from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import plot_fermi_surface
from typing import Union, List, Dict, Any

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

def analyze_band_structure(band_structure_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: analyze_band_structure
    Description: Analyzes the band structure and returns key properties.
    Parameters:
        band_structure_pickle_file_path: str, the path to the pickle file of band structure.
    Returns:
        dict, a dictionary containing key properties of the band structure.
    """
    band_structure = load_from_pickle(band_structure_pickle_file_path)
    
    return {
        'band_gap': band_structure.get_band_gap(),
        'is_metal': band_structure.is_metal(),
        'vbm': band_structure.get_vbm(),
        'cbm': band_structure.get_cbm()
    }

def calculate_density_of_states(dos_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: calculate_density_of_states
    Description: Calculates the density of states and returns key properties.
    Parameters:
        dos_pickle_file_path: str，the path to a pickle file containing the density of states.
    Returns:
        dict, a dictionary containing key properties of the density of states.
    """
    if isinstance(dos_pickle_file_path, str):
        dos = load_from_pickle(dos_pickle_file_path)

    dos_dict = {
        'total_dos': dos.get_densities(),
        'fermi_level': dos.efermi
    } 
    return dos_dict

def plot_fermi_surface(dos_pickle_file_path: str) -> None:
    """
    Name: plot_fermi_surface
    Description: Plots the Fermi surface for the given density of states.
    Parameters:
        dos_pickle_file_path: str, the path to a pickle file containing the density of states.
    Returns:
        None
    """
    if isinstance(dos_pickle_file_path, str):
        dos = load_from_pickle(dos_pickle_file_path)

    fermi_surface_plotter = FermiSurfacePlotter(dos)
    fermi_surface_plotter.get_plot()

