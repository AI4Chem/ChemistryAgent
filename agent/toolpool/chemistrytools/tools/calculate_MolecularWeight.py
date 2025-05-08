from chemistry_tools.pubchem.compound import Compound


def get_compound_MolecularWeight_by_CID(compound_CID: int) -> float:
    """
    Name: get_compound_MolecularWeight_by_CID
    Description: Retrieve the molecular weight (molecular mass) of a chemical compound using its PubChem Compound Identifier (CID).
    Parameters:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    Returns:
        molecular_weight: float, the molecular weight (molecular mass) of the chemical compound.
    """
    compound = Compound.from_cid(compound_CID)
    return compound.molecular_weight


if __name__ == "__main__":
    get_compound_MolecularWeight_by_CID(962)
    get_compound_MolecularWeight_by_CID(5234)
    get_compound_MolecularWeight_by_CID(5988)
