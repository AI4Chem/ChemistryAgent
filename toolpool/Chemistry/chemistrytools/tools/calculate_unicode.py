from chemistry_tools.formulae.unicode import string_to_unicode


def compound_string_to_unicode(compound: str) -> str:
    """
    Name: compound_string_to_unicode
    Description: Convert a chemical formula string to its Unicode representation.
    Parameters:
        compound: str, compound molecular formula in the string format.
    Returns:
        unicode: str, the Unicode representation of the chemical formula.
    """
    unicode = string_to_unicode(compound)
    return unicode


if __name__ == "__main__":
    compound_string_to_unicode("NH4+")
    compound_string_to_unicode("Fe(CN)6+2")
    compound_string_to_unicode("alpha-FeOOH(s)")
