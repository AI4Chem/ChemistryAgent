import pickle
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry
from typing import List, Dict

def save_to_pickle(phase_diagram: PhaseDiagram, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Saves the PhaseDiagram instance to a pickle file.
    Parameters:
        phase_diagram: PhaseDiagram - The PhaseDiagram instance to be saved.
        file_path: str - The path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(phase_diagram, f)

def load_from_pickle(file_path: str) -> PhaseDiagram:
    """
    Name: load_from_pickle
    Description: Loads a PhaseDiagram instance from a pickle file.
    Parameters:
        file_path: str - The path to the pickle file.
    Returns:
        PhaseDiagram - The loaded PhaseDiagram instance.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def create_phase_diagram(entries_pickle_file_path: str, output_phase_diagram_pickle_file_path: str) -> PhaseDiagram:
    """
    Name: create_phase_diagram
    Description: Initializes a PhaseDiagram which loads from a pickle file containing a list of entries then saves it to a pickle file.
    Parameters:
        entries_pickle_file_path: str, the path to the pickle file of computed entries.
        output_phase_diagram_pickle_file_path: str, the path to the pickle file of PhaseDiagram instance to be saved.

    Returns:
        PhaseDiagram, the initialized PhaseDiagram instance.
    """
    entries = load_from_pickle(entries_pickle_file_path)
    phase_diagram = PhaseDiagram(entries)
    save_to_pickle(phase_diagram, output_phase_diagram_pickle_file_path)
    return phase_diagram


def get_decomposition(phase_diagram_pickle_file_path: str, entry_pickle_file_path: str) -> Dict:
    """
    Name: get_decomposition
    Description: Retrieves the decomposition of a given entry in the phase diagram.
    Parameters:
        phase_diagram_pickle_file_path: str, the path to the pickle file of PhaseDiagram instance.
        entry_pickle_file_path: str, the pickle path to entry to get the decomposition for.
    Returns:
        Dict, the decomposition of the entry.
    """
    phase_diagram = load_from_pickle(phase_diagram_pickle_file_path)
    entry = load_from_pickle(entry_pickle_file_path)
    return phase_diagram.get_decomposition(entry)

def get_e_above_hull(phase_diagram_pickle_file_path: str, entry_pickle_file_path: str) -> float:
    """
    Name: get_e_above_hull
    Description: Retrieves the energy above the hull of a given entry in the phase diagram.
    Parameters:
        phase_diagram_pickle_file_path: str, the path to the pickle file of PhaseDiagram instance.
        entry_pickle_file_path: str, the pickle path to entry to get the energy above the hull for.
    Returns:
        float, the energy above the hull of the entry.
    """
    phase_diagram = load_from_pickle(phase_diagram_pickle_file_path)
    entry = load_from_pickle(entry_pickle_file_path)
    return phase_diagram.get_e_above_hull(entry)
