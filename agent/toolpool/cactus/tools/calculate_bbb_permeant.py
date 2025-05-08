"""Tool to calculate the Blood Brain Barrier Permeability of a compound."""

from .adme_pred import ADME
from rdkit import Chem

def CalculateBBBPermeant(compound_SMILES: str) -> bool:
    '''
    Name: CalculateBBBPermeant
    Description: Calculate the Blood Brain Barrier Permeability of the compound.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        judgement: bool, Whether the molecule reaches the blood brain barrier permeant.
    Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    mol = ADME(mol)
    if mol.boiled_egg_bbb():
        return "Yes"
    else:
        return "No"

if __name__ == "__main__":
    CalculateBBBPermeant("CC(C)(C)C(=O)C(N1C=NC=N1)OC2=CC=C(C=C2)Cl")
    CalculateBBBPermeant("CC1=CCC(CC1)C(C)(C)O")
    CalculateBBBPermeant("CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)([O-])OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=C(C7=NC(=CC8=NC(=C(C4=N5)C)C(C8(C)C)CCC(=O)N)C(C7(C)CC(=O)N)CCC(=O)N)C)[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3]")