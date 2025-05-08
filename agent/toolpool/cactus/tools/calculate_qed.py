"""Tool to calculate the QED of a compound."""

from rdkit.Chem import Descriptors, MolFromSmiles

def CalculateQED(compound_SMILES: str) -> float:
    '''
    Name: CalculateQED
    Description: Compute Quantitative Estimate of Druglikeness (QED) of the given molecule.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        QED: float, The QED from 0 (druglike) to 1 (not).
    '''
    return Descriptors.qed(MolFromSmiles(compound_SMILES))

if __name__ == "__main__":
    CalculateQED("CCC(=O)[O-].[Na+]")
    CalculateQED("CC1=C(C=CC=C1C2=CC=CC=C2)COC(=O)C3C(C3(C)C)C=C(C(F)(F)F)Cl")
    CalculateQED("C(C1C2C(C(C(O1)OC3C(OC(C(C3O)O)OC4C(OC(C(C4O)O)OC5C(OC(C(C5O)O)OC6C(OC(C(C6O)O)OC7C(OC(C(C7O)O)OC8C(OC(O2)C(C8O)O)CO)CO)CO)CO)CO)CO)O)O)O")