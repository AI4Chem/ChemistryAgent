from chemlib import Galvanic_Cell
from chemlib import electrolysis

def galvanic_cell_properties(element_electrode1: str, element_electrode2: str) -> dict:
    """
    Name: create_galvanic_cell
    Description: Create a Galvanic (Voltaic) Cell with the given electrodes and return its properties.
    Parameters:
        element_electrode1: str, the elemental composition of one of the electrodes.
        element_electrode2: str, the elemental composition of the other electrode.
    Returns:
        properties: dict, a dictionary of the cell's properties including the cell potential.
    Raises:
        NotImplementedError: If either of the electrodes is invalid or its reduction potential is unknown.
    """
    try:
        cell = Galvanic_Cell(element_electrode1, element_electrode2)
        properties = cell.properties
        return properties
    except NotImplementedError as e:
        return str(e)


def galvanic_cell_potential(element_electrode1: str, element_electrode2: str) -> float:
    """
    Name: create_galvanic_cell
    Description: Create a Galvanic (Voltaic) Cell with the given electrodes and return its potential.
    Parameters:
        element_electrode1: str, the elemental composition of one of the electrodes.
        element_electrode2: str, the elemental composition of the other electrode.
    Returns:
        cell_potential: float, the cell potential of the galvanic cell.
    Raises:
        NotImplementedError: If either of the electrodes is invalid or its reduction potential is unknown.
    """
    try:
        cell = Galvanic_Cell(element_electrode1, element_electrode2)
        cell_potential = cell.cell_potential
        return cell_potential
    except NotImplementedError as e:
        return str(e)


# def galvanic_cell_diagram():
#     return


def perform_electrolysis(element: str, n: int, grams: float, seconds: float) -> dict:
    """
    Name: perform_electrolysis
    Description: Perform electrolysis calculations given an element, the moles of electrons transferred, and two of the following: current (amps), time (seconds), or mass (grams).
    Parameters:
        element: str, the symbol of a chemical element.
        n: int, the moles of electrons transferred.
        grams: float, mass of the meterial.
        seconds: float, time of electrolysis.
    Returns:
        dict: contains the element, moles of electrons transferred, and the calculated values for amps, seconds, and grams.
    Raises:
        TypeError: If not only 2 of the parameters in kwargs are specified.
    """
    result = electrolysis(element, n, grams=grams, seconds=seconds)
    return result

if __name__ == "__main__":
    # Example inputs for testing the function
    perform_electrolysis('Cu', 2, amps=10, seconds=12*60*60)  # Copper metal production
    perform_electrolysis("Ag", 2, amps=2, grams=28.3)  # Electroplating a flute with silver
    perform_electrolysis("Al", 3, grams=805, seconds=24*60*60)  # Aluminum metal production

