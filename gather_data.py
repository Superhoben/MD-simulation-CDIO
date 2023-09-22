from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor


def get_ASE_atoms_from_material_id(material_id: str):
    """Takes a material id as a string and returns the primitive unitcell of the material
    in a format that can be used be ASE simulations"""

    # This uses Gustav Wassbäck's personal API-key to access the database
    with MPRester("Aumz0uNirwQYwJgWgrLVFq3Fr1Z4SfwK") as mpr: 
        some_material = mpr.materials.search(material_ids=[material_id])
        return AseAtomsAdaptor.get_atoms(some_material[-1].structure)


def find_materials_by_elements_and_bandgap(elements: list[str], band_gap: tuple[float, float], fields: list[str]):
    """Takes a list of elements in formula form as strings, band gap range as tuple of two floats (lowest and highest band gap values in eV)
    and a list of proporties as strings. Returns a dictionary with material ids as keys and a doc with the corresponding material properties 
    as values.

    Searches the database for materials containing all the specified elements of parameter elemnts while having a bandgap in the range
    specified by parameter bandgap. Gathers information about only the properties specified by fields."""

    # This uses Gustav Wassbäck's personal API-key to access the database
    if not ("material_id" in fields):
        fields.append("material_id")
    with MPRester("Aumz0uNirwQYwJgWgrLVFq3Fr1Z4SfwK") as mpr: 
        matching_materials = mpr.materials.summary.search(elements=elements, band_gap=band_gap, fields=fields)
        material_dictionary = {}
        for material in matching_materials:
            material_id = str(material.material_id)
            material_dictionary[material_id] = material
        return material_dictionary


# Example to show how it works
materials_dict = find_materials_by_elements_and_bandgap(["Ni", "Sb", "Zr"], (0, 1), ["band_gap"])
print(materials_dict.keys())
print(materials_dict.values())
