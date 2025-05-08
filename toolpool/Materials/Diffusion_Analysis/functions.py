import pickle
from pymatgen.analysis.diffusion.neb import NEBAnalysis
from pymatgen.analysis.diffusion.neb.full_path_mapper import MigrationGraph
from typing import List

def save_to_pickle(data, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Saves the NEBAnalysis or MigrationGraph instance to a pickle file.
    Parameters:
        data: NEBAnalysis or MigrationGraph, the data to be saved.
        file_path: str, the path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str):
    """
    Name: load_from_pickle
    Description: Loads a NEBAnalysis or MigrationGraph instance from a pickle file.
    Parameters:
        file_path: str, the path to the pickle file.
    Returns:
        NEBAnalysis or MigrationGraph, the loaded instance.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def create_neb_analysis(neb_dir: str, output_neb_analysis_pickle_file_path: str) -> NEBAnalysis:
    """
    Name: create_neb_analysis
    Description: Initializes NEBAnalysis from a directory and saves it to a pickle file.
    Parameters:
        neb_dir: str, path to the root directory of the NEB calculation.
        output_neb_analysis_pickle_file_path: str, the path to a pickle file to save the NEBAnalysis instance created.
    Returns:
        NEBAnalysis, the initialized NEBAnalysis instance.
    """

    neb_analysis = NEBAnalysis.from_dir(neb_dir)
    save_to_pickle(neb_analysis, output_neb_analysis_pickle_file_path)
    return neb_analysis

def get_decomposition(neb_analysis_pickle_file_path: str) -> List[float]:
    """
    Name: get_decomposition
    Description: Retrieves the energy decomposition along the NEB path.
    Parameters:
        neb_analysis_pickle_file_path: str, the path to the pickle file containing a NEBAnalysis instance.
    Returns:
        List[float], the energy decomposition along the NEB path.
    """
    return load_from_pickle(neb_analysis_pickle_file_path).get_decomposition()

def get_e_above_hull(neb_analysis_pickle_file_path: str) -> List[float]:
    """
    Name: get_e_above_hull
    Description: Retrieves the energy above the hull for the NEB path.
    Parameters:
        neb_analysis_pickle_file_path: str, the path to the pickle file containing a saved NEBAnalysis instance.
    Returns:
        List[float], the energy above the hull along the NEB path.
    """
    return load_from_pickle(neb_analysis_pickle_file_path).get_e_above_hull()

def create_migration_graph(structure_file_path: str, migrating_ion: str, output_migration_graph_pickle_file_path: str) -> MigrationGraph:
    """
    Name: create_migration_graph
    Description: Initializes MigrationGraph with a structure file and migrating ion then saves it to a pickle file.
    Parameters:
        structure_file_path: str, the path to the structure file.
        migrating_ion: str, the migrating ion species.
        output_migration_graph_pickle_file_path: str, the path to the MigrationGraph file to be saved.
    Returns:
        MigrationGraph, the initialized MigrationGraph instance.
    """
    migration_graph = MigrationGraph.from_structure_file(structure_file_path, migrating_ion)
    save_to_pickle(migration_graph, output_migration_graph_pickle_file_path)
    return migration_graph

def get_migration_paths(migration_graph_pickle_file_path: str) -> List:
    """
    Name: get_migration_paths
    Description: Retrieves possible migration paths from the MigrationGraph.
    Parameters:
        migration_graph_pickle_file_path: str, the path to the MigrationGraph file.
    Returns:
        list, the list of possible migration paths.
    """
    return load_from_pickle(migration_graph_pickle_file_path).get_paths()
