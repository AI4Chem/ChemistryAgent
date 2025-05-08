import pickle
from pymatgen.core.surface import Slab
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.structure import Structure
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

def model_catalyst_surface(structure_pickle_file_path: str, miller_index: List[int], min_slab_size: float, min_vacuum_size: float, output_slab_pickle_file_path: str) -> Slab:
    """
    Name: model_catalyst_surface
    Description: Models the catalyst surface using the given structure and parameters then saves the created slab to a pickle file.
    Parameters:
        structure_pickle_file_path: str, the structure path to its pickle file.
        miller_index: List[int], the Miller index for the surface.
        min_slab_size: float, the minimum size of the slab.
        min_vacuum_size: float, the minimum size of the vacuum.
        output_slab_pickle_file_path: str, path to the pickle file to save the created slab.
    Returns:
        Slab, the created slab.
    """
    if isinstance(structure_pickle_file_path, str):
        structure = load_from_pickle(structure_pickle_file_path)

    slab = Slab(structure, miller_index, min_slab_size, min_vacuum_size)
    save_to_pickle(slab, output_slab_pickle_file_path)
    return slab

def analyze_reaction_pathway(image_pickle_file_path: str, output_pickle_file_path: str) -> NEBAnalysis:
    """
    Name: analyze_reaction_pathway
    Description: Analyzes the reaction pathway using NEBAnalysis and saves is to a pickle file.
    Parameters:
        image_pickle_file_path: str, the path to the pickle file containing the image structure.
        output_pickle_file_path: str, path to the pickle file to save the NEB analysis.
    Returns:
        NEBAnalysis, the NEB analysis result.
    """
    structures = load_from_pickle(image_pickle_file_path)
    neb_analysis = NEBAnalysis.from_endpoints_and_images(structures[0], structures[-1], structures[1:-1])

    save_to_pickle(neb_analysis, output_pickle_file_path)

    return neb_analysis

def calculate_adsorption_energy(slab_pickle_file_path: str, adsorbate: str) -> Dict[str, Any]:
    """
    Name: calculate_adsorption_energy
    Description: Calculates the adsorption energy using AdsorbateSiteFinder.
    Parameters:
        slab_pickle_file_path: str, the path to a pickle file containing the slab.
        adsorbate: str, the adsorbate molecule or atom.
    Returns:
        Dict[str, Any], the adsorption energy calculations.
    """
    if isinstance(slab_pickle_file_path, str):
        slab = load_from_pickle(slab_pickle_file_path)

    adsorbate_structure = Structure.from_file(adsorbate)
    finder = AdsorbateSiteFinder(slab)
    adsorption_sites = finder.find_adsorption_sites()
    adsorption_energies = finder.calculate_adsorption_energies(adsorbate_structure)

    result = {
        'adsorption_sites': adsorption_sites,
        'adsorption_energies': adsorption_energies
    }
    return result
