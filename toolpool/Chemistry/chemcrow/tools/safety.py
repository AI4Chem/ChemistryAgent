import re

import pandas as pd
import pkg_resources
import requests

from ..utils import is_smiles, pubchem_query2smiles, tanimoto


def ExplosiveCheck(compound_CAS: str) -> str:
    '''
    Name: ExplosiveCheck
    Description: Input compound CAS number, returns if molecule is explosive.
    Parameters:
        compound_CAS: str, the CAS number of compound.
    Returns:
        whether molecule is explosive.
    '''
    pubchem_data = {}
    def _fetch_pubchem_data(cas_number):
        """Fetch data from PubChem for a given CAS number, or use cached data if it's already been fetched."""
        if cas_number not in pubchem_data:
            try:
                url1 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/cids/JSON"
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{requests.get(url1).json()['IdentifierList']['CID'][0]}/JSON"
                r = requests.get(url)
                pubchem_data[cas_number] = r.json()
            except:  # noqa: E722
                return "Invalid molecule input, no Pubchem entry."
        return pubchem_data[cas_number]

    def ghs_classification(text):
        """Gives the ghs classification from Pubchem. Give this tool the name or CAS number of one molecule."""
        if is_smiles(text):
            return "Please input a valid CAS number."
        data = _fetch_pubchem_data(text)
        if isinstance(data, str):
            return "Molecule not found in Pubchem."
        try:
            for section in data["Record"]["Section"]:
                if section.get("TOCHeading") == "Chemical Safety":
                    ghs = [
                        markup["Extra"]
                        for markup in section["Information"][0]["Value"][
                            "StringWithMarkup"
                        ][0]["Markup"]
                    ]
                    if ghs:
                        return ghs
        except (StopIteration, KeyError):
            return None

    """Checks if a molecule has an explosive GHS classification using pubchem."""
    # first check if the input is a CAS number
    if is_smiles(compound_CAS):
        return "Please input a valid CAS number."
    cls = ghs_classification(compound_CAS)
    if cls is None:
        return (
            "Explosive Check Error. The molecule may not be assigned a GHS rating. "
        )
    if "Explos" in str(cls) or "explos" in str(cls):
        return "Molecule is explosive"
    else:
        return "Molecule is not known to be explosive"


def SimilarControlChemCheck(compound_SMILES: str) -> str:
    '''
    Name: SimilarControlChemCheck
    Description: Check max similarity between compound and controlled chemicals.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        max similarity with controlled chemicals.
    '''
    def _tanimoto(s1, s2):
        sim = tanimoto(s1, s2)
        if isinstance(sim, float):
            return sim
        return 0.0
    
    data_path = pkg_resources.resource_filename("toolpool", "chemcrow/data/chem_wep_smi.csv")
    cw_df = pd.read_csv(data_path)

    try:
        if not is_smiles(compound_SMILES):
            return "Please input a valid SMILES string."

        max_sim = cw_df["smiles"].apply(lambda x: _tanimoto(compound_SMILES, x)).max()
        if max_sim > 0.35:
            return (
                f"{compound_SMILES} has a high similarity "
                f"({max_sim:.4}) to a known controlled chemical."
            )
        else:
            return (
                f"{compound_SMILES} has a low similarity "
                f"({max_sim:.4}) to a known controlled chemical. "
                "This is substance is safe, you may proceed with the task."
            )
    except:  # noqa: E722
        return "Tool error."


def ControlChemCheck(compound_SMILES: str) -> str:
    '''
    Name: ControlChemCheck
    Description: Input compound in SMILES or CAS number, check if compound is a controlled chemical.
    Parameters:
        compound_SMILES: str, compound in SMILES format.
    Returns:
        whether compound is a controlled chemical.
    '''
    # data_path = pkg_resources.resource_filename("chemcrow", "data/chem_wep_smi.csv")
    data_path = pkg_resources.resource_filename("toolpool", "chemcrow/data/chem_wep_smi.csv")
    cw_df = pd.read_csv(data_path)
    try:
        if is_smiles(compound_SMILES):
            query_esc = re.escape(compound_SMILES)
            found = (
                cw_df["smiles"]
                .astype(str)
                .str.contains(f"^{query_esc}$", regex=True)
                .any()
            )
        else:
            found = (
                cw_df["cas"]
                .astype(str)
                .str.contains(f"^\({compound_SMILES}\)$", regex=True)
                .any()
            )
        if found:
            return (
                f"The molecule {compound_SMILES} appears in a list of "
                "controlled chemicals."
            )
        else:
            # Get smiles of CAS number
            try:
                smi = pubchem_query2smiles(compound_SMILES)
            except ValueError as e:
                return str(e)
            # Check similarity to known controlled chemicals
            return SimilarControlChemCheck(smi)

    except Exception as e:
        return f"Error: {e}"


if __name__ == "__main__":
    ExplosiveCheck("118-96-7")
    ExplosiveCheck("10025-87-3")
    ExplosiveCheck("55-63-0")

    SimilarControlChemCheck("O=P(Cl)(Cl)Cl")
    SimilarControlChemCheck("CC(=O)C")
    SimilarControlChemCheck("O=C1N(C)C(C2=C(N=CN2C)N1C)=O")

    ControlChemCheck("O=P(Cl)(Cl)Cl")
    ControlChemCheck("CC(=O)C")
    ControlChemCheck("O=C1N(C)C(C2=C(N=CN2C)N1C)=O")