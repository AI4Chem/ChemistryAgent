from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..utils import *  # noqa: F403

def MolSimilarity(compound_SMILES_1: str, compound_SMILES_2: str) -> str:
    '''
    Name: MolSimilarity
    Description: Input two molecule SMILES, returns their Tanimoto similarity.
    Parameters:
        compound_SMILES_1: str, compound in SMILES format.
        compound_SMILES_2: str, compound in SMILES format.
    Returns:
        Tanimoto similarity: str.
    '''
    smiles1, smiles2 = compound_SMILES_1, compound_SMILES_2

    similarity = tanimoto(smiles1, smiles2)  # noqa: F405

    if isinstance(similarity, str):
        return similarity

    sim_score = {
        0.9: "very similar",
        0.8: "similar",
        0.7: "somewhat similar",
        0.6: "not very similar",
        0: "not similar",
    }
    if similarity == 1:
        return "Error: Input Molecules Are Identical"
    else:
        val = sim_score[
            max(key for key in sim_score.keys() if key <= round(similarity, 1))
        ]
        message = f"The Tanimoto similarity between {smiles1} and {smiles2} is {round(similarity, 4)},\
        indicating that the two molecules are {val}."
    return message


def SMILES2Weight(compound_SMILES: str) -> str:
    '''
    Name: SMILES2Weight
    Description: Input compound SMILES, returns molecular weight.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        molecular weight: float.
    '''
    mol = Chem.MolFromSmiles(compound_SMILES)
    if mol is None:
        return "Invalid SMILES string"
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    return mol_weight


def FuncGroups(compound_SMILES: str) -> str:
    '''
    Name: FuncGroups
    Description: Input a molecule SMILES or name, return a list of functional groups identified by their common name (in natural language).
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        a list of functional groups identified by their common name (in natural language).
    '''
    # List obtained from https://github.com/rdkit/rdkit/blob/master/Data/FunctionalGroups.txt
    dict_fgs = {
        "furan": "o1cccc1",
        "aldehydes": " [CX3H1](=O)[#6]",
        "esters": " [#6][CX3](=O)[OX2H0][#6]",
        "ketones": " [#6][CX3](=O)[#6]",
        "amides": " C(=O)-N",
        "thiol groups": " [SH]",
        "alcohol groups": " [OH]",
        "methylamide": "*-[N;D2]-[C;D3](=O)-[C;D1;H3]",
        "carboxylic acids": "*-C(=O)[O;D1]",
        "carbonyl methylester": "*-C(=O)[O;D2]-[C;D1;H3]",
        "terminal aldehyde": "*-C(=O)-[C;D1]",
        "amide": "*-C(=O)-[N;D1]",
        "carbonyl methyl": "*-C(=O)-[C;D1;H3]",
        "isocyanate": "*-[N;D2]=[C;D2]=[O;D1]",
        "isothiocyanate": "*-[N;D2]=[C;D2]=[S;D1]",
        "nitro": "*-[N;D3](=[O;D1])[O;D1]",
        "nitroso": "*-[N;R0]=[O;D1]",
        "oximes": "*=[N;R0]-[O;D1]",
        "Imines": "*-[N;R0]=[C;D1;H2]",
        "terminal azo": "*-[N;D2]=[N;D2]-[C;D1;H3]",
        "hydrazines": "*-[N;D2]=[N;D1]",
        "diazo": "*-[N;D2]#[N;D1]",
        "cyano": "*-[C;D2]#[N;D1]",
        "primary sulfonamide": "*-[S;D4](=[O;D1])(=[O;D1])-[N;D1]",
        "methyl sulfonamide": "*-[N;D2]-[S;D4](=[O;D1])(=[O;D1])-[C;D1;H3]",
        "sulfonic acid": "*-[S;D4](=O)(=O)-[O;D1]",
        "methyl ester sulfonyl": "*-[S;D4](=O)(=O)-[O;D2]-[C;D1;H3]",
        "methyl sulfonyl": "*-[S;D4](=O)(=O)-[C;D1;H3]",
        "sulfonyl chloride": "*-[S;D4](=O)(=O)-[Cl]",
        "methyl sulfinyl": "*-[S;D3](=O)-[C;D1]",
        "methyl thio": "*-[S;D2]-[C;D1;H3]",
        "thiols": "*-[S;D1]",
        "thio carbonyls": "*=[S;D1]",
        "halogens": "*-[#9,#17,#35,#53]",
        "t-butyl": "*-[C;D4]([C;D1])([C;D1])-[C;D1]",
        "tri fluoromethyl": "*-[C;D4](F)(F)F",
        "acetylenes": "*-[C;D2]#[C;D1;H]",
        "cyclopropyl": "*-[C;D3]1-[C;D2]-[C;D2]1",
        "ethoxy": "*-[O;D2]-[C;D2]-[C;D1;H3]",
        "methoxy": "*-[O;D2]-[C;D1;H3]",
        "side-chain hydroxyls": "*-[O;D1]",
        # "side-chain aldehydes": "*=[O;D1]",
        "primary amines": "*-[N;D1]",
        "nitriles": "*#[N;D1]",
    }
    def _is_fg_in_mol(mol, fg):
        fgmol = Chem.MolFromSmarts(fg)
        mol = Chem.MolFromSmiles(mol.strip())
        return len(Chem.Mol.GetSubstructMatches(mol, fgmol, uniquify=True)) > 0
    try:
        fgs_in_molec = [
            name
            for name, fg in dict_fgs.items()
            if _is_fg_in_mol(compound_SMILES, fg)
        ]
        if len(fgs_in_molec) > 1:
            return f"This molecule contains {', '.join(fgs_in_molec[:-1])}, and {fgs_in_molec[-1]}."
        else:
            return f"This molecule contains {fgs_in_molec[0]}."
    except:  # noqa: E722
        return "Wrong argument. Please input a valid molecular SMILES."


if __name__ == "__main__":
    MolSimilarity("O=C1N(C)C(C2=C(N=CN2C)N1C)=O.CC(C)c1ccccc1")
    MolSimilarity("O=C1N(C)C(C2=C(N=CN2C)N1C)=O.O=C1N(C)C(C2=C(N=CN2C)N1CCC)=O")
    MolSimilarity("O=C1N(C)C(C2=C(N=CN2C)N1C)=O")
    # MolSimilarity("O=C1N(C)C(C2=C(N=CN2C)N1C)=O.4-(4-hydroxyphenyl)butan-2-one")

    SMILES2Weight("O=C1N(C)C(C2=C(N=CN2C)N1C)=O")
    SMILES2Weight("CC(C)c1ccccc1")
    SMILES2Weight("O=C1N(C)C(C2=C(N=CN2C)N1C)=Ox")

    FuncGroups("O=C1N(C)C(C2=C(N=CN2C)N1C)=O")
    FuncGroups("CCCCCCCCC[NH+]1C[C@@H]([C@H]([C@@H]([C@H]1CO)O)O)O")
    FuncGroups("4-(4-hydroxyphenyl)butan-2-one")