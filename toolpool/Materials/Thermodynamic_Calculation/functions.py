import pickle
from pymatgen.core import Composition
from pymatgen.core.units import FloatWithUnit
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
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

def calculate_reaction_energy(reactant_entries_pickle_file_path: str, product_entries_pickle_file_path: str) -> dict:
    """
    Name: calculate_reaction_energy
    Description: Calculates the reaction energy for a given set of reactants and products and saves to a pickle file.
    Parameters:
        reactant_entries_pickle_file_path: str, the path to pickle files containing chemical formula entries representing the reactants.
        product_entries_pickle_file_path: str, the path to pickle files containing chemical formula entries representing the products.
    Returns:
        dict
    """
    
    reactant_entries = load_from_pickle(reactant_entries_pickle_file_path)
    product_entries = load_from_pickle(product_entries_pickle_file_path)

    reaction = ComputedReaction(reactant_entries, product_entries)
    energy = FloatWithUnit(reaction.calculated_reaction_energy, "eV atom^-1").to("kJ mol^-1")

    result = {
        "reaction": str(reaction),
        "reaction_energy": str(energy)
    }
    return result

def get_most_stable_entry(entries_pickle_file_path: str, composition_str: str, output_entry_pickle_file_path: str) -> None:
    """
    Name: get_most_stable_entry
    Description: Retrieves the most stable entry for a given chemical formula from a list of computed entries and saves to a pickle file.
    Parameters:
        entries_pickle_file_path: str, the path to the input pickle file containing the list of computed entries.
        composition_str: str, the chemical formula for which to find the most stable entry.
        output_entry_pickle_file_path: str, the path to the output pickle file where the most stable entry will be saved.
    Returns:
        None
    """
    entries: List[ComputedEntry] = load_from_pickle(entries_pickle_file_path)
    relevant_entries = [
        entry
        for entry in entries
        if entry.composition.reduced_formula == Composition(composition_str).reduced_formula
    ]
    relevant_entries = sorted(relevant_entries, key=lambda e: e.energy_per_atom)
    most_stable_entry = relevant_entries[0]
    save_to_pickle(most_stable_entry, output_entry_pickle_file_path)

def analyze_phase_stability(entries_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: analyze_phase_stability
    Description: Analyzes the phase stability according to the entries and returns key properties.
    Parameters:
        entries_pickle_file_path: str, the path to the input pickle file containing the list of computed entries.
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

def create_pourbaix_diagram(entries_pickle_file_path: str, output_pourbaix_diagram_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: create_pourbaix_diagram
    Description: Creates a Pourbaix diagram then saves it to a pickle file and returns key properties.
    Parameters:
        entries_pickle_file_path: str, the path to the input pickle file containing the list of computed entries.
        output_pourbaix_diagram_pickle_file_path: str, path to the pickle file to save the Pourbaix diagram object.
    Returns:
        dict, a dictionary containing key properties of the Pourbaix diagram.
    """
    entries = load_from_pickle(entries_pickle_file_path)
    pourbaix_diagram = PourbaixDiagram(entries)
    save_to_pickle(pourbaix_diagram, output_pourbaix_diagram_pickle_file_path)
    
    return {
        'stable_entries': pourbaix_diagram.stable_entries,
        'unstable_entries': pourbaix_diagram.unstable_entries,
        'pourbaix_diagram': pourbaix_diagram
    }
