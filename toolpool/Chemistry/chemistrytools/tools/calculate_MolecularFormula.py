from chemistry_tools.pubchem.compound import Compound


def convert_compound_CID_to_Molecular_Formula(compound_CID: int) -> str:
    """
    Name: convert_compound_CID_to_Molecular_Formula
    Description: Retrieve the molecular formula of a chemical compound using its PubChem Compound Identifier (CID).
    Parameters:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    Returns:
        compound: str, the molecular formula of the chemical compound.
    """
    compound = Compound.from_cid(compound_CID)
    return compound.molecular_formula


if __name__ == "__main__":
    convert_compound_CID_to_Molecular_Formula(962)
    convert_compound_CID_to_Molecular_Formula(5234)
    convert_compound_CID_to_Molecular_Formula(5988)
