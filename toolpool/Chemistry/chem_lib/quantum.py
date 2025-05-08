from chemlib import Element, rydberg, Wave
from chemlib import energy_of_hydrogen_orbital as _energy_of_hydrogen_orbital

def electromagnetic_wave_by_wavelength(wavelength: float) -> dict:
    """
    Name: electromagnetic_wave_by_wavelength
    Description: Get the properties of the wave based on the wavelength.
    Parameters:
        wavelength: float, wavelength of the wave.
    Returns:
        properties: dict, properties of the wave (wavelength, frequency, energy).
    """
    return Wave(wavelength=wavelength).properties


def electromagnetic_wave_by_frequency(frequency: float) -> dict:
    """
    Name: electromagnetic_wave_by_wavelength
    Description: Get the properties of the wave based on the frequency.
    Parameters:
        frequency: float, frequency of the wave.
    Returns:
        properties: dict, properties of the wave (wavelength, frequency, energy).
    """
    return Wave(frequency=frequency).properties


def electromagnetic_wave_by_energy(energy: float) -> dict:
    """
    Name: electromagnetic_wave_by_wavelength
    Description: Get the properties of the wave based on the energy.
    Parameters:
        energy: float, energy of the wave.
    Returns:
        properties: dict, properties of the wave (wavelength, frequency, energy).
    """
    return Wave(energy=energy).properties


def electromagnetic_wave_by_rydberg_equation(element: str, n1: int, n2: int) -> dict:
    """
    Name: conduct_rydberg
    Description: Get the properties of the wave based on the Rydberg equation for a given element and orbitals.
    Parameters:
        element: str, the symbol of the element. For example, "H"
        n1: int, initial orbital number.
        n2: int, final orbital number (must be greater than n1).
    Returns:
        properties: dict, properties of the wave (wavelength, frequency, energy) as calculated using Rydberg equation.
    Raises:
        ValueError: If n2 is not greater than n1.
    """
    element = Element(element)
    wave_properties  = rydberg(element, n1, n2)
    return wave_properties


def energy_of_hydrogen_orbital(n: int) -> float:
    """
    Name: calculate_hydrogen_orbital_energy
    Description: Calculate the energy of an electron in the nth orbital of the Hydrogen atom in Joules.
    Parameters:
        n: int, the orbital number (1 for the first orbital, 2 for the second, etc.).
    Returns:
        energy: float, Energy of the electron in Joules.
    """
    energy = _energy_of_hydrogen_orbital(n)
    return energy

if __name__ == "__main__":
    # Example inputs for testing the function
    energy_of_hydrogen_orbital(3)
    energy_of_hydrogen_orbital(2)
    energy_of_hydrogen_orbital(10)

    electromagnetic_wave_by_rydberg_equation("C", 2, 5)

