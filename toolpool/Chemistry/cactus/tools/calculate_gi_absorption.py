"""Tool to calculate the GI Absorption of the compound."""

from .adme_pred import ADME
from rdkit import Chem

def CalculateGIAbsorption(compound_SMILES: str) -> str:
    '''
    Name: CalculateGIAbsorption
    Description: Calculate the GastroIntestinal Absorption (GI Absorption) of the compound.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        GI_Absorption: str, Return whether the GastroIntestinal Absorption (GI Absorption) of a compound is high or low.
    Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    mol = ADME(mol)
    if mol.boiled_egg_hia():
        return "High"
    else:
        return "Low"

if __name__ == "__main__":
    CalculateGIAbsorption("C#C")
    CalculateGIAbsorption("CCCCCCCCCCCCCCCC[N+]1=CC=CC=C1.[Cl-]")
    CalculateGIAbsorption("C1C(C(C(C(C1N)OC2C(C(C(C(O2)CN)O)O)O)O)OC3C(C(C(C(O3)CO)O)N)O)N")