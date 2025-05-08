import pickle
import pymatgen.core as mg

def save_to_pickle(data: any, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Save data to a pickle file.
    Parameters:
        data: any, the data to be saved.
        file_path: str, the path to the output pickle file where the data will be saved.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> any:
    """
    Name: load_from_pickle
    Description: Load data from a pickle file.
    Parameters:
        file_path: str, the path to the input pickle file from which the data will be loaded.
    Returns:
        data: any, the data loaded from the pickle file.
    """
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def get_atomic_mass(element_str: str) -> float:
    """
    Name: get_atomic_mass
    Description: Retrieves the atomic mass of the given element.
    Parameters:
        element_str: str, the symbol of the given element.
    Returns:
        float, the atomic mass of the given element.
    """
    element_symbol = element_str
    element = mg.Element(element_symbol)
    atomic_mass = element.atomic_mass
    return atomic_mass

def get_melting_point(element_str: str) -> float:
    """
    Name: get_melting_point
    Description: Retrieves the melting point of the given element.
    Parameters:
        element_str: str, the symbol of the given element.
    Returns:
        float, the melting point of the given element.
    """
    element_symbol = load_from_pickle(element_str)
    element = mg.Element(element_symbol)
    melting_point = element.melting_point
    return melting_point

def get_ionic_radii(element_str: str) -> dict:
    """
    Name: get_ionic_radii
    Description: Retrieves the ionic radii of the given element.
    Parameters:
        element_str: str, the symbol of the given element.
    Returns:
        dict, all ionic radii of the given element as a dict of {oxidation state: ionic radii}.
    """
    element_symbol = element_str
    element = mg.Element(element_symbol)
    ionic_radii = element.ionic_radii
    return ionic_radii

def get_atomic_mass_in_unit(element_str: str, unit: str) -> float:
    """
    Name: get_atomic_mass_in_unit
    Description: Retrieves the atomic mass of a specified element in a given unit.
    Parameters:
        element_str: str, the symbol of the given element.
        unit: str, the desired unit for the atomic mass (e.g., 'kg', 'g', etc.).
    Returns:
        float, the atomic mass in the specified unit.
    """
    element_symbol = element_str
    element = mg.Element(element_symbol)
    atomic_mass_in_unit = element.atomic_mass.to(unit)
    return atomic_mass_in_unit

def get_species_properties(element_str: str, oxidation_state: int) -> dict:
    """
    Name: get_species_properties
    Description: Retrieves properties of a species with a given oxidation state.
    Parameters:
        element_str: str, the symbol of the given element.
        oxidation_state: int, the oxidation state of the species.
    Returns:
        dict, the species properties.
    """
    element_symbol = element_str
    species = mg.Species(element_symbol, oxidation_state)
    species_properties = {
        "atomic_mass": species.atomic_mass,
        "ionic_radius": species.ionic_radius
    }
    return species_properties

def get_composition_properties(composition_str: str, element_str: str) -> dict:
    """
    Name: get_composition_properties
    Description: Retrieves properties of a specified element in a composition given its formula.
    Parameters:
        composition_str: str, the chemical composition.
        element_str: str, the symbol of the given element.
    Returns:
        dict, the composition properties.
    """
    compound = composition_str
    element = element_str
    comp = mg.Composition(compound)
    composition_properties = {
        "weight": comp.weight,
        "amount_of_element": comp[element],
        "atomic_fraction": comp.get_atomic_fraction(element),
        "weight_fraction": comp.get_wt_fraction(element)
    }
    return composition_properties