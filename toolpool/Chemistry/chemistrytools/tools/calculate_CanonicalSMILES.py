from chemistry_tools.pubchem.compound import Compound


def convert_compound_CID_to_SMILES(compound_CID: int) -> str:
    """
    Name: convert_compound_CID_to_SMILES
    Description: Retrieve the Canonical SMILES representation of a chemical compound using its PubChem Compound Identifier (CID).
    Parameters:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    Returns:
        compound_SMILES: str, the Canonical SMILES representation of the chemical compound.
    """
    compound = Compound.from_cid(compound_CID)
    return compound.smiles


if __name__ == "__main__":
    convert_compound_CID_to_SMILES(962)
    convert_compound_CID_to_SMILES(5234)
    convert_compound_CID_to_SMILES(5988)
