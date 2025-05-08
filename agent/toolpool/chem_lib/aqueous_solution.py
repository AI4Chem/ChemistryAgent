import chemlib
from chemlib import Solution

def acidity_calculation(pH: float) -> dict:
    """
    Name: acidity_calculation
    Description: For any inputted pH, pOH, [H+], or [OH-], finds the corresponding values.
    Parameters:
        pH: The value of the chosen input (pH=, pOH=, H=, or OH=).
    Returns:
        acidity: dict, PH value.
    """
    if type(pH) is not float:
        return "Input Error"
    return chemlib.pH(pH=pH)

def make_solution_by_grams_per_liter(compound: str, grams: float, liters: float) -> float:
    """
    Name: make_solution_by_grams_per_liter
    Description: Given grams of solute per liter. Return the molarity of the solution.
    Parameters:
        compound: str, the formula of the solute compound.
        grams: float, how many grams of solute.
        liters: float, how many liters of solution.
    Returns:
        molarity: float, the molarity of the solution.
    """
    # Calculate molarity from grams per liter
    # Instantiate Solution object
    solution = Solution.by_grams_per_liters(compound, grams,liters)
    return solution.molarity


if __name__ == "__main__":
   make_solution_by_grams_per_liter("NaCl", 10, 1)
   make_solution_by_grams_per_liter("KCl", 10, 10)
   make_solution_by_grams_per_liter("MgCl2", 10, 10)