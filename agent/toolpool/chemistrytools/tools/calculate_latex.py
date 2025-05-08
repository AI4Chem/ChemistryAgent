from chemistry_tools.formulae.latex import string_to_latex


def compound_string_to_latex(compound: str) -> str:
    """
    Name: compound_string_to_latex
    Description: Convert a chemical formula string to its LaTeX representation.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        latex: str, the LaTeX representation of the chemical formula.
    """
    latex = string_to_latex(compound)
    return latex


if __name__ == "__main__":
    compound_string_to_latex("NH4+")
    compound_string_to_latex("Fe(CN)6+2")
    compound_string_to_latex("alpha-FeOOH(s)")
