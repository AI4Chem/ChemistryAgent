"""Tool to calculate the LogP of a compound."""

from rdkit import Chem
from rdkit.Chem import Descriptors

def CalculateLogP(compound_SMILES: str) -> float:
    '''
    Name: CalculateLogP
    Description: Calculate the log of the partition coefficient (LogP) of a compound.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        LogP: float, The log of the partition coefficient (LogP) of a compound.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    return Descriptors.MolLogP(mol)

if __name__ == "__main__":
    CalculateLogP("CC(=O)CCC(=O)O")
    CalculateLogP("C1=CC(=CC=C1[N+](=O)[O-])Cl")
    CalculateLogP("CC(=O)N1CCN(CC1)C2=CC=C(C=C2)OCC3COC(O3)(CN4C=CN=C4)C5=C(C=C(C=C5)Cl)Cl")