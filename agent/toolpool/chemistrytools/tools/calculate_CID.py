from chemistry_tools.pubchem.lookup import get_compounds


def get_compound_CID(compound: str) -> int:
    """
    Name: get_compound_CID
    Description: Retrieve the PubChem Compound Identifier (CID) for a given chemical name.
    Parameters:
        compound: str, compound name in the string format.
    Returns:
        compound_CID: int, the PubChem Compound Identifier (CID) of the chemical compound.
    """
    CID = get_compounds(compound)[0].cid
    return CID


if __name__ == "__main__":
    get_compound_CID("water")
    get_compound_CID("salt")
    get_compound_CID("sugar")
