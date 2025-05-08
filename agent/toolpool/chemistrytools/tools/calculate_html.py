from chemistry_tools.formulae.html import string_to_html

def compound_string_to_html(compound: str) -> str:
    """
    Name: compound_string_to_html
    Description: Convert a chemical formula string to its HTML representation.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        html: str, the HTML representation of the chemical formula.
    """
    html = string_to_html(compound)
    return html


if __name__ == "__main__":
    compound_string_to_html("NH4+")
    compound_string_to_html("Fe(CN)6+2")
    compound_string_to_html("alpha-FeOOH(s)")
