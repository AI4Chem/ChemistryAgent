from chemistry_tools.formulae.formula import Formula


def calculate_compound_molar_mass(compound: str) -> float:
    """
    Name: calculate_compound_molar_mass
    Description: Calculate the molar mass of a given chemical formula.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        mass: float, the molar mass of the chemical formula.
    """
    formula = Formula.from_string(compound)
    return formula.mass


if __name__ == "__main__":
    calculate_compound_molar_mass("H20")
    calculate_compound_molar_mass("H2")
    calculate_compound_molar_mass("O2")
