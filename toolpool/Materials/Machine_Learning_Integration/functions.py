import pickle
import numpy as np
from matminer.featurizers.site import CrystalNNFingerprint
from pymatgen.core.structure import Structure
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from typing import List, Dict, Any, Union

def save_to_pickle(obj: Any, file_path: str) -> None:
    """
    Saves an object to a pickle file.
    
    Parameters:
        obj: Any - The object to be saved.
        file_path: str - The path to the pickle file.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(obj, f)

def load_from_pickle(file_path: str) -> Any:
    """
    Loads an object from a pickle file.
    
    Parameters:
        file_path: str - The path to the pickle file.
    
    Returns:
        Any - The loaded object.
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def generate_features(structure_pickle_file_path_list: List[str], output_feature_pickle_file_path: str) -> np.ndarray:
    """
    Name: generate_features
    Description: Generates features for a list of structures using CrystalNNFingerprint then saves it to a pickle file.
    Parameters:
        structure_pickle_file_path_list: List[Union[Structure, str]], the list of paths to structure pickle files.
        output_feature_pickle_file_path: str, path to the pickle file to save the generated features.
    Returns:
        np.ndarray, the array of generated features.
    """
    fingerprint = CrystalNNFingerprint.from_preset('cn')
    features = []

    for structure in structure_pickle_file_path_list:
        if isinstance(structure, str):
            structure = load_from_pickle(structure)
        feature = fingerprint.featurize(structure)
        features.append(feature)

    features = np.array(features)

    save_to_pickle(features, output_feature_pickle_file_path)

    return features

def train_model(feature_pickle_file_path: str, targets: List[float], test_size: float, random_state: int, output_model_pickle_file_path: str) -> Dict[str, Any]:
    """
    Name: train_model
    Description: Trains a RandomForest model on the provided features and targets then saves it to a pickle file.
    Parameters:
        feature_pickle_file_path: str, the array of features or path to the pickle file.
        targets: List[float], the list of target values.
        test_size: float, the proportion of the dataset to include in the test split.
        random_state: int, the random state for train-test splitting.
        output_model_pickle_file_path: str, path to the pickle file to save the trained model.
    Returns:
        dict, a dictionary containing the trained model, test predictions, and mean squared error.
    """
    if isinstance(feature_pickle_file_path, str):
        features = load_from_pickle(feature_pickle_file_path)

    X_train, X_test, y_train, y_test = train_test_split(features, targets, test_size=test_size, random_state=random_state)
    
    model = RandomForestRegressor()
    model.fit(X_train, y_train)
    
    predictions = model.predict(X_test)
    mse = mean_squared_error(y_test, predictions)

    save_to_pickle(model, output_model_pickle_file_path)

    return {
        'model': model,
        'predictions': predictions,
        'mse': mse
    }

def predict_with_model(model_pickle_file_path: str, feature_pickle_file_path: str) -> np.ndarray:
    """
    Name: predict_with_model
    Description: Makes predictions using a trained model and provided features.
    Parameters:
        model_pickle_file_path: str, path to the pickle file of the trained model.
        feature_pickle_file_path: str, the array of features or path to the pickle file.
    Returns:
        np.ndarray, the array of predictions.
    """
    if isinstance(model_pickle_file_path, str):
        model = load_from_pickle(model_pickle_file_path)
    if isinstance(feature_pickle_file_path, str):
        features = load_from_pickle(feature_pickle_file_path)

    predictions = model.predict(features)
    return predictions
