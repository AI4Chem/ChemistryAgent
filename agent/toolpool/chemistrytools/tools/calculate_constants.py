from chemistry_tools.constants import (
    avogadro_number,
    boltzmann_constant,
    electron_radius,
    faraday_constant,
    molar_gas_constant,
    neutron_mass,
    plancks_constant,
    speed_of_light,
    vacuum_permittivity,
)


def get_constants(name: str) -> tuple:
    """
    Name: get_constants
    Description: Retrieve the value, unit, and symbol of a physical constant by its name.
    Parameters:
        name: str, the name of the physical constant (e.g., 'avogadro_number', 'boltzmann_constant').
    Returns:
        tuple: (name, value, unit, symbol), where
            name: str, the name of the constant.
            value: float, the numerical value of the constant.
            unit: str, the unit of the constant.
            symbol: str, the symbol or abbreviation of the constant.
    """
    constant_dict = {
        "avogadro_number": avogadro_number,
        "boltzmann_constant": boltzmann_constant,
        "electron_radius": electron_radius,
        "faraday_constant": faraday_constant,
        "molar_gas_constant": molar_gas_constant,
        "neutron_mass": neutron_mass,
        "plancks_constant": plancks_constant,
        "speed_of_light": speed_of_light,
        "vacuum_permittivity": vacuum_permittivity,
    }
    constant = constant_dict[name]
    value = constant.value
    unit = constant.unit
    symbol = constant.symbol
    return [name, value, str(unit), symbol]


if __name__ == "__main__":
    get_constants("vacuum_permittivity")
    get_constants("avogadro_number")
    get_constants("boltzmann_constant")
