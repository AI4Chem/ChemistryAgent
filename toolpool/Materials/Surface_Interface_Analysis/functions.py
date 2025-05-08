import pickle
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.structure import Structure
from typing import List, Union, Dict, Any

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

def generate_surface(structure_pickle_file_path: str, miller_index: List[int], min_slab_size: float, min_vacuum_size: float) -> Dict[str, Any]:
    """
    Name: generate_surface
    Description: Generates a surface slab for a given structure.
    Parameters:
        structure_pickle_file_path: Union[Structure, str], the structure object or path to the pickle file.
        miller_index: List[int], the Miller index of the surface.
        min_slab_size: float, minimum size of the slab.
        min_vacuum_size: float, minimum size of the vacuum layer.
    Returns:
        dict, a dictionary containing the slab object and its properties.
    """
    structure = load_from_pickle(structure_pickle_file_path)
    
    slab_generator = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size)
    slab = slab_generator.get_slabs()[0]  # Get the first slab for simplicity
    return {
        'slab': slab,
        'miller_index': slab.miller_index,
        'surface_area': slab.surface_area
    }

def model_interfacial_reactions(entries_pickle_file_path_1:  str, entries_pickle_file_path_2: str, open_elem: str, temperature: float) -> Dict[str, Any]:
    """
    Name: model_interfacial_reactions
    Description: Models interfacial reactions between two sets of entries.
    Parameters:
        entries_pickle_file_path_1: str, path to the pickle file containing the first set of computed entries.
        entries_pickle_file_path_2: str, path to the pickle file containing the second set of computed entries.
        open_elem: str, the open element in the reaction.
        temperature: float, the temperature of the reaction.
    Returns:
        dict, a dictionary containing the interfacial reactivity object and its properties.
    """
    entries1 = load_from_pickle(entries_pickle_file_path_1)
    entries2 = load_from_pickle(entries_pickle_file_path_2)
    interfacial_reactivity = InterfacialReactivity(entries1, entries2, open_elem, temperature)
    return {
        'interfacial_reactivity': interfacial_reactivity,
        'reaction_entries': interfacial_reactivity.get_kinks(),
        'reaction_energies': interfacial_reactivity.get_energies()
    }

def perform_adsorption_study(structure_pickle_file_path: str, adsorbate: str) -> Dict[str, Any]:
    """
    Name: perform_adsorption_study
    Description: Performs an adsorption study on a given structure.
    Parameters:
        structure_pickle_file_path: str, path to the pickle file containing the structure object.
        adsorbate: str, the adsorbate molecule or atom.
    Returns:
        dict, a dictionary containing the adsorption sites and their properties.
    """
    structure = load_from_pickle(structure_pickle_file_path)
    
    adsorbate_site_finder = AdsorbateSiteFinder(structure)
    adsorption_sites = adsorbate_site_finder.find_adsorption_sites(adsorbate=adsorbate)
    
    return {
        'adsorption_sites': adsorption_sites,
        'number_of_sites': len(adsorption_sites['all_sites'])
    }
