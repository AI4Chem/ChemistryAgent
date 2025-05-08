from .safety import ControlChemCheck
from ..utils import (
    is_smiles,
    pubchem_query2smiles,
    query2cas,
)


def Query2CAS(compound_SMILES: str) -> str:
    '''
    Name: Query2CAS
    Description: Input molecule (name or SMILES), returns CAS number.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        compound_CAS: str, the CAS number of compound.
    '''
    url_cid = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/cids/JSON"
    )
    url_data = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON"
    )
    try:
        # if query is smiles
        smiles = None
        if is_smiles(compound_SMILES):
            smiles = compound_SMILES
        try:
            cas = query2cas(compound_SMILES, url_cid, url_data)
        except ValueError as e:
            return str(e)
        if smiles is None:
            try:
                smiles = pubchem_query2smiles(cas, None)
            except ValueError as e:
                return str(e)
        # check if mol is controlled
        msg = ControlChemCheck(smiles)
        if "high similarity" in msg or "appears" in msg:
            return f"CAS number {cas}found, but " + msg
        return cas
    except ValueError:
        return "CAS number not found"

if __name__ == "__main__":
    Query2CAS("4-(4-hydroxyphenyl)butan-2-one") # 5471-51-2
    Query2CAS("O=C1N(C)C(C2=C(N=CN2C)N1C)=O") # 58-08-2
    Query2CAS("nomol") # invalid