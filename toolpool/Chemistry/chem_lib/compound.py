from typing import Literal

from chemlib import Compound

def calculate_element_frequencies_in_compound(compound: str) -> dict:
    """
    Name: calculate_element_frequencies_in_compound
    Description: Calculate the frequencies of the constituent elements in the compound.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        element_frequencies: dict, the frequencies of the constituent elements in the compound.
    """
    return Compound(compound).occurences


def calculate_compound_molar_mass(compound: str) -> float:
    """
    Name: calculate_compound_molar_mass
    Description: Calculate the molar mass in (g/mol) of the compound.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        molar_mass: float, the molar mass in (g/mol) of the compound.
    """
    return Compound(compound).molar_mass()


def calculate_compound_percentage_composition_by_mass(compound: str, element: str) -> float:
    """
    Name: calculate_compound_percentage_composition_by_mass
    Description: Get the percentage composition by mass of a certain element of the compound.
    Parameters:
        compound: str, compound molecular formula in the string format.
        element: str, the constituent element of which the user wants to get percentage composition.
    Returns:
        element_percentage_compostion: float, the percentage composition by mass of the element in the compound.
    """
    return Compound(compound).percentage_by_mass(element)


def convert_compound_stoichiometry_amount(compound: str, unit: Literal["grams", "moles", "molecules"], amount: float) -> dict:
    """
    Name: convert_compound_stoichiometry_amount
    Description: Get stoichiometric amounts of the compound given one measurement.
    Parameters:
        compound: str, compound molecular formula in the string format.
        unit: Literal["grams", "moles", "molecules"], given stoichiometry amount.
        amount: float, the amount of the chosen compound.
    Returns:
        compound_stoichiometry_amount: dict, the percentage composition by mass of the element in the compound.
    """
    return eval(f"Compound('{compound}').get_amounts({unit}={amount})")


if __name__ == "__main__":
    calculate_element_frequencies_in_compound()