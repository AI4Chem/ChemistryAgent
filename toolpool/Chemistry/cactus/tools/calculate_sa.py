"""Tool to calculate the SA of a compound."""

from rdkit.Chem import MolFromSmiles
from . import sascorer

def CalculateSA(compound_SMILES: str) -> float:
    '''
    Name: CalculateSA
    Description: Used to compute the synthetic accessibility (SA) of the given molecule in the given SMILES string.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        SA: float, The SA between 1 (easy) and 10 (hard).
    '''
    return sascorer.calculateScore(MolFromSmiles(compound_SMILES))

if __name__ == "__main__":
    CalculateSA("CC1=C(C(=O)CC1OC(=O)C2C(C2(C)C)C=C(C)C)CC=C")
    CalculateSA("C(=O)(O)[O-].[Na+]")
    CalculateSA("C1=CC(=CC=C1[N+](=O)[O-])Cl")
