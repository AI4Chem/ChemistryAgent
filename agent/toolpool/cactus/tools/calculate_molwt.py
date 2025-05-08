"""Tool to calculate the MolWt of a compound."""

from rdkit.Chem import Descriptors, MolFromSmiles

def CalculateMolWt(compound_SMILES: str) -> float:
    '''
    Name: CalculateMolWt
    Description: Calculate the molecular weight (MolWt) of the given compound in SMILES string. Units in Dalton.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        molecular_weight: float, The exact molecular weight in daltons.
    '''
    mol = MolFromSmiles(compound_SMILES)
    return Descriptors.ExactMolWt(mol)

if __name__ == "__main__":
    CalculateMolWt("C1=CC=C(C=C1)N")
    CalculateMolWt("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34")
    CalculateMolWt("CC(=O)N1CCN(CC1)C2=CC=C(C=C2)OCC3COC(O3)(CN4C=CN=C4)C5=C(C=C(C=C5)Cl)Cl")