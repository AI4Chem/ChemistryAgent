from pymatgen.core import Element

def create_element(element_str: str) -> Element:
    """
    Name: create_element
    Description: Initializes an Element with a given symbol.
    Parameters:
        element_str: str, the symbol of the element (e.g., 'Si' for Silicon).    
    Returns:
        Element, the initialized Element instance.
    """
    return Element(element_str)

def get_atomic_mass(element_str: str) -> float:
    """
    Name: get_atomic_mass
    Description: Retrieves the atomic mass of the element.
    Parameters:
        element_str: str, the symbol of the element (e.g., 'Si' for Silicon).   
    Returns:
        float, the atomic mass of the element.
    """
    return Element(element_str).atomic_mass

def get_atomic_radius(element_str: str) -> float:
    """
    Name: get_atomic_radius
    Description: Retrieves the atomic radius of the element.
    Parameters:
        element_str: str, the symbol of the element (e.g., 'Si' for Silicon).   
    Returns:
        float, the atomic radius of the element.
    """
    return Element(element_str).atomic_radius

def get_electron_affinity(element_str: str) -> float:
    """
    Name: get_electron_affinity
    Description: Retrieves the electron affinity of the element.
    Parameters:
        element_str: str, the symbol of the element (e.g., 'Si' for Silicon).   
    Returns:
        float, the electron affinity of the element.
    """
    return Element(element_str).electron_affinity
