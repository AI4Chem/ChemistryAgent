import pickle
from typing import Any, List, Dict
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.apps.battery.analyzer import BatteryAnalyzer
from pymatgen.core.composition import Composition
from pymatgen.ext.matproj import MPRester

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

def get_all_entries(chemical_system: List[str], output_entries_pickle_file_path: str) -> None:
    """
    Name: get_all_entries
    Description: Retrieves all entries for the specified chemical system from the Materials Project database and saves to a pickle file.
    Parameters:
        chemical_system: List[str], a list of elements defining the chemical system (e.g., ["Ca", "C", "O"]).
        output_entries_pickle_file_path: List[str], the path to the output pickle file where the entries will be saved.
    Returns:
        None
    """
    
    with MPRester() as mpr:
        entries = mpr.get_entries_in_chemsys(chemical_system)
    save_to_pickle(entries, output_entries_pickle_file_path)

def balance_reaction(reactant_entries_pickle_file_path: str, product_entries_pickle_file_path: str, output_reaction_pickle_file_path: str) -> str:
    """
    Name: balance_reaction
    Description: Balances a chemical reaction with reactant and product entries pickle file then saves it to a pickle file.
    Parameters:
        reactant_entries_pickle_file_path: str, the path to pickle files containing chemical formula entries representing the reactants.
        product_entries_pickle_file_path: str, the path to pickle files containing chemical formula entries representing the products.
        output_reaction_pickle_file_path: str, path to the pickle file containing the reaction to be saved.
    Returns:
        str, the balanced reaction as a string.
    """
    reactants = load_from_pickle(reactant_entries_pickle_file_path)
    products = load_from_pickle(product_entries_pickle_file_path)

    reaction = Reaction(reactants, products)
    save_to_pickle(reaction, output_reaction_pickle_file_path)
    return str(reaction)

def analyze_battery(structure_pickle_file_path: str) -> dict:
    """
    Name: analyze_battery
    Description: Performs battery-related calculations.
    Parameters:
        structure_pickle_file_path: str, path to the pickle file containing a saved Composition structure instance.
    Returns:
        dict, a dictionary of battery analysis results.
    """
    structure = load_from_pickle(structure_pickle_file_path)
    
    analyzer = BatteryAnalyzer(structure)
    return analyzer.analyze()
