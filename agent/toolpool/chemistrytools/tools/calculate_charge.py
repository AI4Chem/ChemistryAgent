from chemistry_tools.pubchem.compound import Compound


def get_compound_charge_by_CID(compound_CID: int) -> int:
    """
    Name: get_compound_charge_by_CID
    Description: Retrieve the charge of a chemical compound using its PubChem Compound Identifier (CID).
    Parameters:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    Returns:
        charge: int, the charge of the chemical compound.
    """
    compound = Compound.from_cid(compound_CID)
    return compound.charge


if __name__ == "__main__":
    get_compound_charge_by_CID(962)
    get_compound_charge_by_CID(5234)
    get_compound_charge_by_CID(70678590)
