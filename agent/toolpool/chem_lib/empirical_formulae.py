from chemlib import empirical_formula_by_percent_comp as efbpc
from chemlib.thermochemistry import combustion_analysis



def get_empirical_formula_by_percent_composition(percent_composition: dict) -> str:
    """
    Name: get_empirical_formula_by_percent_composition
    Description: Calculate the empirical formula of a compound based on the given percentage compositions of elements.
    Parameters:
        percent_composition: dict, keyword arguments representing the percentage compositions of elements in the compound (e.g., {'C'=80.6, 'H'=19.4}).
    Returns:
        compound_unicode: str, the empirical formula of the compound in unicode format.
    Raises:
        ValueError: If the sum of the percentages is not equal to 100.
    """
    if sum(percent_composition.values()) != 100:
        raise ValueError("The sum of the percentage compositions must be equal to 100.")

    empirical_formula = efbpc(**percent_composition).formula
    return empirical_formula


def analyze_combustion(CO2: float, H2O: float) -> str:
    """
    Name: analyze_combustion
    Description: Get the empirical formula of a hydrocarbon given the grams of CO2 and grams of H2O formed from its combustion.
    Parameters:
        CO2: float, the grams of carbon dioxide formed as a result of the combustion of the hydrocarbon.
        H2O: float, the grams of water formed as a result of the combustion of the hydrocarbon.
    Returns:
        compound_unicode: str, the empirical formula of the hydrocarbon.
    """
    return combustion_analysis(CO2, H2O)


if __name__ == "__main__":
    # Example inputs for testing the function
    get_empirical_formula_by_percent_composition({"C":80.6, "H":19.4})  # Example: Compound with 80.6% C and 19.4% H
    get_empirical_formula_by_percent_composition({"N":30.4, "O":69.6})  # Example: Compound with 30.4% N and 69.6% O
    get_empirical_formula_by_percent_composition({"S":50.0, "O":50.0})  # Example: Compound with 50% S and 50% O


