"""Tool to calculate the TPSA of a compound."""

from rdkit.Chem import MolFromSmiles, rdMolDescriptors

def CalculateTPSA(compound_SMILES: str) -> float:
    '''
    Name: CalculateTPSA
    Description: Compute the Topological polar surface area (TPSA) of the given molecule.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        TPSA: float, The TPSA in angstroms^2.
    '''
    return rdMolDescriptors.CalcTPSA(MolFromSmiles(compound_SMILES))

if __name__ == "__main__":
    CalculateTPSA("C1CC(NC1)C(=O)O")
    CalculateTPSA("C1=CC(=CC=C1C(=O)O)C(=O)O")
    CalculateTPSA("CCCCCCCCCCCCCCCC(=O)OCC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C")