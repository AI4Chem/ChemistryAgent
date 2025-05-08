import pickle
from pymatgen.core.surface import SlabGenerator, Structure, Lattice
from pymatgen.analysis.wulff import WulffShape

def save_to_pickle(data, file_path: str) -> None:
    """
    Name: save_to_pickle
    Description: Saves the data to a pickle file.
    Parameters:
        data: any, The data to be saved.
        file_path: str, The path to the pickle file.
    Returns:
        None
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path: str) -> any:
    """
    Name: load_from_pickle
    Description: Loads data from a pickle file.
    Parameters:
        file_path: str, The path to the pickle file.
    Returns:
        data: any, The loaded data.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def get_slabs_info(structure_pickle_file_path: str):
    """
    Name: def get_slabs_info(structure_file_path: str):
    Description: Get slabs information using SlabGenerator with the provided parameters.
    Parameters:
        structure_pickle_file_path: str, path to the pickle file containing a Structure instance.
    Returns:
        slabs: list, a list of generated slab structures.
    """
    slab_structure = load_from_pickle(structure_pickle_file_path)
    slab_generator = SlabGenerator(slab_structure)
    return slab_generator.get_slabs()

def analyze_wulff_shape(lattice_pickle_file_path: str, miller_energies: dict, output_wulff_shape_pickle_file_path: str):
    """
    Name: analyze_wulff_shape
    Description: Analyzes the equilibrium crystal shape using WulffShape and saves WulffShape to a pickle file.
    Parameters:
        lattice_pickle_file_path: str, path to the pickle file containing a Lattice instance.
        miller_energies: dict, dictionary mapping Miller indices to surface energies.
        output_wulff_shape_pickle_file_path: str, path to pickle file containing wulff shape to be saved.
    Returns:
        wulff_shape: WulffShape, The WulffShape object.
    """
    wulff_lattice = load_from_pickle(lattice_pickle_file_path)
    wulff_shape = WulffShape(wulff_lattice, miller_energies)
    save_to_pickle(wulff_shape, output_wulff_shape_pickle_file_path)
    return wulff_shape
