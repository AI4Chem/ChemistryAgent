from chemistry_tools.pubchem.compound import Compound


def convert_compound_CID_to_IUPAC(compound_CID: int) -> str:
    """
    Name: convert_compound_CID_to_IUPAC
    Description: Retrieve the IUPAC name of a chemical compound using its PubChem Compound Identifier (CID).
    Parameters:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    Returns:
        iupac_name: str, the IUPAC name of the chemical compound.
    """
    compound = Compound.from_cid(compound_CID)
    return compound.iupac_name


if __name__ == "__main__":
    convert_compound_CID_to_IUPAC(962)
    convert_compound_CID_to_IUPAC(5234)
    convert_compound_CID_to_IUPAC(5988)
