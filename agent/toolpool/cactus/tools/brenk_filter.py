"""Tool for calculating if a compound passes the Brenk Filter."""

from .adme_pred import ADME
from rdkit import Chem

def BrenkFilter(compound_SMILES: str) -> bool:
    '''
    Name: BrenkFilter
    Description: Calculate whether a molecule triggers the Brenk Filter.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        judgement: bool, Boolean of whether the molecule triggers the Brenk filter.
    Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    mol = ADME(mol)
    return mol.brenk()

if __name__ == "__main__":
    BrenkFilter("CCOP(=S)(OCC)OC1=CC=C(C=C1)[N+](=O)[O-]")
    BrenkFilter("C1CCOC1")
    BrenkFilter("C1=CC=C2C(=C1)C(=CC=C2S(=O)(=O)[O-])N=NC3=C4C=CC(=CC4=CC(=C3O)S(=O)(=O)[O-])S(=O)(=O)[O-].[Na+].[Na+].[Na+]")