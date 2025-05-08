from chemistry_tools.elements._elements import ELEMENTS


def get_element_information(element: str) -> dict:
    """
    Name: get_element_information
    Description: Calculate information about the given element, including its atomic number, symbol, name, atomic mass, description, electron configuration, and electron configuration dictionary.
    Parameters:
        element: a str, the atomic symbol of an element.
    Returns:
        number: int, the atomic number of the element.
        symbol: str, the symbol of the element.
        name: str, the name of the element.
        mass: float, the atomic mass of the element.
        description: str, a brief description of the element.
        eleconfig: str, the electron configuration of the element.
        eleconfig_dict: dict, the electron configuration in dictionary form.
    """
    ele = ELEMENTS[element]
    number = ele.number
    symbol = ele.symbol
    name = ele.name
    mass = ele.mass
    description = ele.description
    eleconfig = ele.eleconfig
    eleconfig_dict = str(ele.eleconfig_dict)

    return {
        "number": number,
        "symbol": symbol,
        "name": name,
        "mass": mass,
        "description": description,
        "eleconfig": eleconfig,
        "eleconfig_dict": eleconfig_dict
    }


if __name__ == "__main__":
    get_element_information("C")
    get_element_information("H")
    get_element_information("O")
