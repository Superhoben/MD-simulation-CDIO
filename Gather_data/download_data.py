"""Material gathering module.

This module contains functions to gather data from the next-gen materials project database using
the so called new API (as opposed to the legacy API).
"""

from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.trajectory import Trajectory
from ase.build import make_supercell
from numpy import cbrt
from math import floor
import os.path


def make_traj_from_material_id(material_id: str, target_number_of_atoms=8000):
    """Take a material id and save a corresponding atoms object in
    a .traj file which can be found in the folder Trajectory files.

    Args:
        material_id (str): The material id must exist in the next-gen material
            projects database and should be in the form 'mp-1234' but with
            another number.
        traget_number_of_atoms (int): The number of atoms of the material desired
            in the system. It is a maximum rather than a minimum.

    Returns:
        none
    """
    # This uses Gustav Wassbäck's personal API-key to access the database
    with MPRester("Aumz0uNirwQYwJgWgrLVFq3Fr1Z4SfwK") as mpr:
        some_material = mpr.materials.search(material_ids=[material_id])
        primitive_cell = AseAtomsAdaptor.get_atoms(some_material[-1].structure)
        number_atoms_primitive = primitive_cell.get_number_of_atoms()
        n = floor(cbrt(target_number_of_atoms/number_atoms_primitive))
        M = [[n, 0, 0], [0, n, 0], [0, 0, n]]
        atoms = make_supercell(primitive_cell, M)
        size_descripition = str(atoms.get_number_of_atoms())+'_atoms_of_'
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
        location_and_name = path_to_traj_folder + size_descripition + material_id + '.traj'
        traj = Trajectory(location_and_name, "w")
        traj.write(atoms)


def get_ASE_atoms_from_material_id(material_id: str):
    """Take a material id and return the primitive unitcell of
    the material in a format that can be used in ASE simulations.

    Args:
        material_id (str): The material id must exist in the next-gen material
            projects database and should be in the form 'mp-1234' but with
            another number.

    Returns:
        atoms: An ASE atoms object, for more information see
            https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    """
    # This uses Gustav Wassbäck's personal API-key to access the database
    with MPRester("Aumz0uNirwQYwJgWgrLVFq3Fr1Z4SfwK") as mpr:
        some_material = mpr.materials.search(material_ids=[material_id])
        return AseAtomsAdaptor.get_atoms(some_material[-1].structure)


def find_materials_by_elements_and_bandgap(elements: list[str], band_gap: tuple[float, float], fields: list[str]):
    """Search for materials containing at least every material mentioned in 'elements' that has a band gap in the
    range specified by 'band_gap' and for these materials gather only the properties specified by 'fields'.

    Args:
        elements (list[str]): Elements are stated in chemical short notation as found in the periodic table. An
            example argument could be ['Cu', 'Fe']. Only elements 1 to 83 and 89 to 94 can be used as of 26/9-2023.
        band_gap (tuple[float, float]): Filter for materials with a bandgap in the specified range, where the first
            float is the lowest allowed bandgap value and the second float is the highest in eV. Due to no bandgap
            being less than 0 and the highest known bandgap being less than 15eV a requirement for the argument is
            0 ≤ First float ≤ Second float ≤ 15. An example argument when looking for conductors could be (0, 0.3),
            when looking for semi conductors (0.6, 2) and when looking for insulators (5, 15).
        fields (list[str]): Specifies which properties that will be gathered for the found materials. The properties
            can be read about at https://docs.materialsproject.org/methodology/materials-methodology and every valid
            field name can be seen in the example argument below:
            ['builder_meta', 'nsites', 'elements', 'nelements', 'composition', 'composition_reduced', 'formula_pretty',
            'formula_anonymous', 'chemsys', 'volume', 'density', 'density_atomic', 'symmetry', 'property_name',
            'material_id', 'deprecated', 'deprecation_reasons', 'last_updated', 'origins', 'warnings', 'structure',
            'task_ids', 'uncorrected_energy_per_atom', 'energy_per_atom', 'formation_energy_per_atom', 'energy_above_hull',
            'is_stable', 'equilibrium_reaction_energy_per_atom', 'decomposes_to', 'xas', 'grain_boundaries', 'band_gap',
            'cbm', 'vbm', 'efermi', 'is_gap_direct', 'is_metal', 'es_source_calc_id', 'bandstructure', 'dos', 'dos_energy_up',
            'dos_energy_down', 'is_magnetic', 'ordering', 'total_magnetization', 'total_magnetization_normalized_vol',
            'total_magnetization_normalized_formula_units', 'num_magnetic_sites', 'num_unique_magnetic_sites',
            'types_of_magnetic_species', 'k_voigt', 'k_reuss', 'k_vrh', 'g_voigt', 'g_reuss', 'g_vrh', 'universal_anisotropy',
            'homogeneous_poisson', 'e_total', 'e_ionic', 'e_electronic', 'n', 'e_ij_max', 'weighted_surface_energy_EV_PER_ANG2',
            'weighted_surface_energy', 'weighted_work_function', 'surface_anisotropy', 'shape_factor', 'has_reconstructed',
            'possible_species', 'has_props', 'theoretical', 'database_IDs']

    Returns:
        dict[str, MPDataDoc]: A dictionary with material id's a keys and MPDataDocs object containing all properties for the
            corresponding material that was asked for by the fields argument.
    """
    if not ("material_id" in fields):
        fields.append("material_id")
    # This uses Gustav Wassbäck's personal API-key to access the database
    with MPRester("Aumz0uNirwQYwJgWgrLVFq3Fr1Z4SfwK") as mpr:
        matching_materials = mpr.materials.summary.search(elements=elements, band_gap=band_gap, fields=fields)
        material_dictionary = {}
        for material in matching_materials:
            material_id = str(material.material_id)
            material_dictionary[material_id] = material
        return material_dictionary


# Example to show how it works
if __name__ == "__main__":
    materials_dict = find_materials_by_elements_and_bandgap(["Ni", "Sb", "Zr"], (0, 1), ["band_gap"])
    print(materials_dict.keys())
    print(materials_dict.values())
    make_traj_from_material_id('mp-24')
