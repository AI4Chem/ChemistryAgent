from chemlib import Element

def get_element_properties(element:str) -> dict:
    """
    Name: get_element_properties
    Description: Get all the properties of the given element, including its atomic number, name, symbol, atomic mass, neutrons, protons, electrons, period, group, phase, whether it's radioactive, whether it's natural, whether it's a metal, nonmetal or metalloid, type, atomic radius, electronegativity, first ionization energy, density, melting point, boiling point, isotopes, discoverer, year of discovery, specific heat, shells, valence, electron configuration, and mass number.
    Parameters:
        element: str, the atomic symbol of an element.
    Returns:
            AtomicNumber: float, the atomic number of the element.
            Element: str, the name of the element.
            Symbol: str, the symbol of the element.
            AtomicMass: float, the atomic mass of the element.
            Neutrons: float, the number of neutrons in the element.
            Protons: float, the number of protons in the element.
            Electrons: float, the number of electrons in the element.
            Period: float, the period in the periodic table where the element is found.
            Group: float, the group in the periodic table where the element is found.
            Phase: str, the phase of the element at room temperature.
            Radioactive: bool, whether the element is radioactive.
            Natural: bool, whether the element occurs naturally.
            Metal: bool, whether the element is a metal.
            Nonmetal: bool, whether the element is a nonmetal.
            Metalloid: bool, whether the element is a metalloid.
            Type: str, the type of the element (e.g., Noble Gas).
            AtomicRadius: str, the atomic radius of the element.
            Electronegativity: float, the electronegativity of the element.
            FirstIonization: str, the first ionization energy of the element.
            Density: str, the density of the element.
            MeltingPoint: str, the melting point of the element.
            BoilingPoint: str, the boiling point of the element.
            Isotopes: float, the number of isotopes of the element.
            Discoverer: str, the discoverer of the element.
            Year: str, the year the element was discovered.
            SpecificHeat: str, the specific heat of the element.
            Shells: float, the number of electron shells in the element.
            Valence: float, the valence electrons in the element.
            Config: str, the electron configuration of the element.
            MassNumber: float, the mass number of the element.
    """
    ele = Element(element)
    AtomicNumber = ele.AtomicNumber
    ElementName = ele.Element
    Symbol = ele.Symbol
    AtomicMass = ele.AtomicMass
    Neutrons = ele.Neutrons
    Protons = ele.Protons
    Electrons = ele.Electrons
    Period = ele.Period
    Group = ele.Group
    Phase = ele.Phase
    Radioactive = ele.Radioactive
    Natural = ele.Natural
    Metal = ele.Metal
    Nonmetal = ele.Nonmetal
    Metalloid = ele.Metalloid
    Type = ele.Type
    AtomicRadius = ele.AtomicRadius
    Electronegativity = ele.Electronegativity
    FirstIonization = ele.FirstIonization
    Density = ele.Density
    MeltingPoint = ele.MeltingPoint
    BoilingPoint = ele.BoilingPoint
    Isotopes = ele.Isotopes
    Discoverer = ele.Discoverer
    Year = ele.Year
    SpecificHeat = ele.SpecificHeat
    Shells = ele.Shells
    Valence = ele.Valence
    Config = ele.Config
    MassNumber = ele.MassNumber
    return {
        'AtomicNumber': AtomicNumber,
        'Element': ElementName,
        'Symbol': Symbol,
        'AtomicMass': AtomicMass,
        'Neutrons': Neutrons,
        'Protons': Protons, 
        'Electrons': Electrons, 
        'Period': Period,
        'Group': Group,
        'Phase': Phase,
        'Radioactive': Radioactive,
        'Natural': Natural,
        'Metal': Metal,
        'Nonmetal': Nonmetal,
        'Metalloid': Metalloid,
        'Type': Type,
        'AtomicRadius': AtomicRadius,
        'Electronegativity': Electronegativity,
        'FirstIonization': FirstIonization,
        'Density': Density,
        'MeltingPoint': MeltingPoint,
        'BoilingPoint': BoilingPoint,
        'Isotopes': Isotopes,
        'Discoverer': Discoverer,
        'Year': Year,
        'SpecificHeat': SpecificHeat,
        'Shells': Shells,
        'Valence': Valence,
        'Config': Config,
        'MassNumber': MassNumber
    }

if __name__ == "__main__":
    get_element_properties("Xe")  # Example: Xenon
    get_element_properties("C")   # Example: Carbon
    get_element_properties("O")   # Example: Oxygen