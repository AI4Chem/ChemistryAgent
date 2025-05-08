import molbloom

from ..utils import is_multiple_smiles, split_smiles


def PatentCheck(compound_SMILES: str) -> str:
    """
    Name: PatentCheck
    Description: Input compound in SMILES format, returns if molecule is patented. You may also input several compounds, separated by a period.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        whether molecule is patented.
    """
    if is_multiple_smiles(compound_SMILES):
        smiles_list = split_smiles(compound_SMILES)
    else:
        smiles_list = [compound_SMILES]
    try:
        output_dict = {}
        for smi in smiles_list:
            r = molbloom.buy(smi, canonicalize=True, catalog="surechembl")
            if r:
                output_dict[smi] = "Patented"
            else:
                output_dict[smi] = "Novel"
        return str(output_dict)
    except:  # noqa: E722
        return "Invalid SMILES string"

if __name__ == "__main__":
    PatentCheck("O=C1N(C)C(C2=C(N=CN2C)N1C)=O") # Patented
    PatentCheck("CCCCCCCCC[NH+]1C[C@@H]([C@H]([C@@H]([C@H]1CO)O)O)O") # Novel
    PatentCheck("4-(4-hydroxyphenyl)butan-2-one") # Invalid SMILES string (It's a molecule with iupac name.)
