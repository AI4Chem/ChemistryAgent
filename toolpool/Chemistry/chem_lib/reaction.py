from typing import Literal

from chemlib import Compound, Reaction, Combustion

def judge_the_balance_of_reaction(reaction_formula: str) -> bool:
    """
    Name: judge_the_balance_of_reaction
    Description: Judge whether the reaction formula is balanced.
    Parameters:
        reaction_formula: str, a string representing the reaction formula.
    Returns:
        Is_Balanced: bool, whether the formula is balanced.
    """
    reaction = Reaction.by_formula(reaction_formula)
    Is_Balanced = reaction.is_balanced
    return Is_Balanced


def reactant_formulas_of_reaction(reaction_formula: str) -> list:
    """
    Name: reactant_formulas_of_reaction
    Description: Get reactant formulas of the reaction.
    Parameters:
        reaction_formula: str, a string representing the reaction formula.
    Returns:
        Reactant_Formulas: list, the reactant formulas.
    """
    reaction = Reaction.by_formula(reaction_formula)
    Reactant_Formulas = reaction.reactant_formulas
    return Reactant_Formulas


def product_formulas_of_reaction(reaction_formula: str) -> list:
    """
    Name: product_formulas_of_reaction
    Description: Get product formulas of the reaction.
    Parameters:
        reaction_formula: str, a string representing the reaction formula.
    Returns:
        Product_Formulas: list, the product formulas.
    """
    reaction = Reaction.by_formula(reaction_formula)
    Product_Formulas = reaction.product_formulas
    return Product_Formulas


def combustion_reactions(compound: str) -> tuple[str, bool, list, list]:
    """
    Name: combustion_reactions
    Description: Create and balance a combustion reaction for a given compound.
    Parameters:
        compound: str, the molecular formula of the compound to combust.
    Returns:
        Formula: str, the balanced combustion reaction formula.
        Is_Balanced: bool, whether the formula is balanced.
    """
    compound = Compound(compound)
    combustion_reaction = Combustion(compound)
    Formula = combustion_reaction.formula
    Is_Balanced = combustion_reaction.is_balanced
    Reactant_Formulas = combustion_reaction.reactant_formulas
    Product_Formulas = combustion_reaction.product_formulas
    return Formula, Is_Balanced # , Reactant_Formulas, Product_Formulas


def balance_the_reaction(reaction_formula: str) -> str:
    """
    Name: balance_the_reaction
    Description: Balance the chemical reaction.
    Parameters:
        reaction_formula: str, a string representing the reaction formula.
    Returns:
        reaction_formula: str, the balanced reaction formula.
        Is_Balanced: bool, whether the formula is balanced.
        Reactant_Formulas: list, the reactant formulas.
        Product_Formulas: list, the product formulas.
    """
    reaction = Reaction.by_formula(reaction_formula)
    reaction.balance()
    Formula = reaction.formula
    return Formula


def reaction_stoichiometry_amounts(reaction_formula: str, compound_number: int, unit: Literal["grams", "moles", "molecules"], amount: float) -> list:
    """
    Name: reaction_stoichiometry_amounts
    Description: Get stoichiometric amounts of all compounds in the reaction given the amount of one compound.
    Parameters:
        reaction_formula: str, the reaction formula string.
        compound_number: int, the chosen compound in the reaction by order of appearance (left to right).
        unit: Literal["grams", "moles", "molecules"], given stoichiometry amount.
        amount: float, the amount of the chosen compound.
    Returns:
        List of dicts: each dict contains the amounts of each compound in the reaction (Compound, Grams, Moles, Molecules).
    Raises:
        ValueError: if the compound_number is less than 1 or greater than the number of compounds in the reaction.
        ValueError: if more than one argument is given under kwargs.
    """
    reaction = Reaction.by_formula(reaction_formula)
    return eval(f"reaction.get_amounts({compound_number}, {unit}={amount})")


def limiting_reagent_of_reaction(reaction_formula: str, mode: list[Literal["grams", "moles", "molecules"]], rectant_amount_list: list) -> list:
    """
    Name: find_limiting_reagent
    Description: Determine the limiting reagent in a chemical reaction.
    Parameters:
        reaction_formula: str, the reaction formula string.
        mode: list, the units of each amount in args. Default is grams, can also be moles or molecules.
        rectant_amount_list: list, the amounts of each reactant to use in the chemical reaction.
    Returns:
        Limiting_Reagent: str, the formula of the limiting reagent.
    Raises:
        TypeError: If the number of args doesn’t match the number of reactants in the reaction.
        ValueError: If the mode is not grams, moles, or molecules.
    """
    reaction = Reaction.by_formula(reaction_formula)
    limiting_reagent = reaction.limiting_reagent(*rectant_amount_list, mode=mode)
    return limiting_reagent.formula

if __name__ == "__main__":
    # combustion_reactions("CH4")  # Methane
    # combustion_reactions("C2H6")  # Ethane
    # combustion_reactions("C3H8")  # Propane

    # reaction_stoichiometry_amounts("N2O5 + H2O --> HNO3", 1, "grams", 5)  # Given 5 grams of N2O5
    # reaction_stoichiometry_amounts("N2O5 + H2O --> HNO3", 3, "moles", 3.5)  # Given 3.5 moles of HNO3
    # reaction_stoichiometry_amounts("H2 + O2 --> H2O", 1, "grams", 2)  # Given 2 grams of H2

    limiting_reagent_of_reaction("N2O5 + H2O --> HNO3", 'moles', [50, 80])  # Given 50 grams of N2O5 and 80 grams of H2O
    limiting_reagent_of_reaction("N2O5 + H2O --> HNO3", 'moles', [3, 1])  # Given 3 moles of N2O5 and 1 mole of H2O
    limiting_reagent_of_reaction("H2 + O2 --> H2O", 'grams', [2, 1])  # Given 2 grams of H2 and 1 gram of O2