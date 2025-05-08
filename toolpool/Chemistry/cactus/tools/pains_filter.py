"""Tool for calculating if a molecule passes the PAINS Filter."""

from .adme_pred import ADME
from rdkit import Chem

def PainsFilter(compound_SMILES: str) -> bool:
    '''
    Name: PainsFilter
    Description: Calculate whether a molecule triggers the Pains Filter (pan-assay interference compounds).
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        judgement: bool, Boolean of whether the molecule triggers the PAINS filter.
    Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    mol = ADME(mol)
    return mol.pains()

if __name__ == "__main__":
    PainsFilter("CC1=CC=CC=C1")
    PainsFilter("CC(=CCC1=C(C=CC2=C1OC(=O)C=C2)OC)C")
    PainsFilter("CC(=O)N1CCN(CC1)C2=CC=C(C=C2)OCC3COC(O3)(CN4C=CN=C4)C5=C(C=C(C=C5)Cl)Cl")