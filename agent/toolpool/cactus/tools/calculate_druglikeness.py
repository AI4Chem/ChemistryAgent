"""Tool to calculate if a compound passes the Lipinski rule of five."""

from .adme_pred import ADME
from rdkit import Chem


def CalculateDruglikeness(compound_SMILES: str) -> str | list | bool:
    '''
    Name: CalculateDruglikeness
    Description: Calculate the druglikeness of the compound with regards to Lipinski's Rule of Five.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        judgement: The druglikeness of the compound with regards to Lipinski's Rule of Five.
    Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    mol = ADME(mol)
    return mol.druglikeness_lipinski(verbose=True)

if __name__ == "__main__":
    CalculateDruglikeness("CSSC")
    CalculateDruglikeness("C(SC#N)SC#N")
    CalculateDruglikeness("C1C(C(C(C(C1N)OC2C(C(C(C(O2)CN)O)O)O)O)OC3C(C(C(C(O3)CO)O)N)O)N")
