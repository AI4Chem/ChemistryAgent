from typing import List, Dict
from pymatgen.ext.matproj import MPRester

def get_element_properties(element_str: str) -> Dict[str, float]:
    """
    Name: get_element_properties
    Description: Retrieves properties of a specified element from the Materials Project database.
    Parameters:
        element_str: str, the symbol of the element (e.g., 'Au' for Gold).
    Returns:
        dict, a dictionary containing properties such as atomic mass, density, and electron affinity.
    Other Information: Requires an API key from Materials Project.
    """
    # Replace 'YOUR_API_KEY' with your actual Materials Project API key
    API_KEY = 'YOUR_API_KEY'
    properties = {}

    with MPRester(API_KEY) as mpr:
        data = mpr.query({"symbol": element_str}, ["material_id", "atomic_mass", "density", "electron_affinity"])
        if data:
            properties = {
                "material_id": data[0]["material_id"],
                "atomic_mass": data[0]["atomic_mass"],
                "density": data[0]["density"],
                "electron_affinity": data[0]["electron_affinity"]
            }
    return properties

def calculate_average(number_list: List[float]) -> float:
    """
    Name: calculate_average
    Description: Calculates the average of a list of numbers.
    Parameters:
        number_list: List[float], a list of floating-point numbers.
    Returns:
        float, the average of the numbers.
    Other Information: Returns 0.0 for an empty list.
    """
    if not number_list:
        return 0.0
    return sum(number_list) / len(number_list)
