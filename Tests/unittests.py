import os.path, sys, unittest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.build import bulk, molecule
from ase import Atoms, Atom
#from asap3 import EMT
from ase.calculators.emt import EMT
from tkinter import Tk
from Simulation.simple_simulation import run_simple_md_simulation
from Simulation.lattice_constant import optimize_scaling
from Simulation.calc_properties import approx_lattice_constant, calc_temp, calc_pressure, calc_mean_square_displacement, lindemann_criterion, self_diffusion_coefficent
from Simulation.calc_bulk_properties import calc_bulk_modulus, calculate_cohesive_energy
from Simulation.run_md_simulation import run_single_md_simulation
from Gather_data.download_data import get_ASE_atoms_from_material_id
from User_interface.user_interface import initiate_gui
from User_API_key import start_program
import numpy as np
from ase import units


class UnitTests(unittest.TestCase):
    """Definitions of unittests.

    Different tests can be added here from different files. To add a new test,
    import the functions from the file that is to be tested and then write a 
    new test function below the existing ones.

    See example tests below.
    """
    def test_calc_temp(self):
        atoms = FaceCenteredCubic(symbol="Cu", size=(5, 5, 5), pbc=True)
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)
        self.assertTrue(270 < calc_temp(atoms) < 330)

    def test_lattice_constant_gradient_descent(self):
        """Uses a material with known lattice constant, disorts it and see if the lattice method is able
        to find the original lattice constant which also should be the optimal one."""
        #atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms = FaceCenteredCubic(symbol="Cu", size=(2, 2, 2), pbc=True)
        atoms.set_cell(atoms.cell*0.5, scale_atoms=True)  # Rescale unit cell so atoms are now in suboptimal positions
        atoms.calc = EMT()
        scaling = optimize_scaling(atoms, {'optimal_scaling': [], 'iterations_to_find_scaling': []})
        # print(scaling)
        self.assertTrue((0.98 < scaling*0.5) and (scaling*0.5 < 1.02))
        

    # More test are needed for this function
    def test_calc_pressure_no_field(self):
        #atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms = FaceCenteredCubic(symbol="Cu", size=(2, 2, 2), pbc=True)
        atoms.calc = EMT()  # If atoms are only, Ni, Cu, Pd, Ag, Pt or Au no input parameters are nessecary
        # When atoms are still at optimal positions no pressure should occur, +-0.1 GPa tolerance is given
        self.assertTrue((-0.1 < calc_pressure(atoms)) and (calc_pressure(atoms) < 0.1))
        
    def test_bulk_modulus(self):
        # lattice_constant = 4
        # choosing more specific silver to test the bulk modulus
        #atoms = get_ASE_atoms_from_material_id('mp-124')
        atoms = FaceCenteredCubic(symbol="Ag", size=(2, 2, 2), pbc=True)
        atoms.calc = EMT()
        # Another way to create our object "bulk" this time
        #ag_bulk = bulk(name= "Ag",crystalstructure= "fcc", a=4.09)
        #ag_bulk.calc = EMT
        calc_bulk_modulus(atoms)
        # From material website https://next-gen.materialsproject.org/materials/mp-124?chemsys=Ag#elastic_constants
        # we have the bulk modulus for Ag to be 88 GPa
        # According to Rickard, it is absolutely no problem to get 100 GPa of bulk modulus since we are first taking the
        # second derivative approximation and secondly we are using EMT calculator.
        self.assertTrue((86 < calc_bulk_modulus(atoms)) and (calc_bulk_modulus(atoms) < 105))

    def test_cohesive_energy(self):
        # Cohesive energy per atom (eV/atom) values from Charles Kittle book "Introduction to Solid State Physics" page 50
        # Ag cohesive energy = 2.95
        # Au cohesive energy = 3.81
        # Ni cohesive energy = 4.44
        # Cu cohesive energy = 3.49
        # Create multiple atoms object
        atom_structure_Ag = Atoms("Ag")
        atom_structure_Ag.calc = EMT()
        atom_structure_Au = Atoms("Au")
        atom_structure_Au.calc = EMT()
        atom_structure_Ni = Atoms("Ni")
        atom_structure_Ni.calc = EMT()
        atom_structure_Cu = Atoms("Cu")
        atom_structure_Cu.calc = EMT()
        # From 
        # you can find every single cohesive energy per atom for each elements in 
        self.assertTrue((2.9 < calculate_cohesive_energy(atom_structure_Ag)) and
                        (calculate_cohesive_energy(atom_structure_Ag) < 3) and 
                        (3.76 < calculate_cohesive_energy(atom_structure_Au)) and 
                        (calculate_cohesive_energy(atom_structure_Au) < 3.86) and 
                        (4.39 < calculate_cohesive_energy(atom_structure_Ni)) and
                        (calculate_cohesive_energy(atom_structure_Ni) < 4.49) and 
                        (3.44 < calculate_cohesive_energy(atom_structure_Cu)) and
                        (calculate_cohesive_energy(atom_structure_Cu) < 3.54))

    def test_cohesive_energy2(self):
        # Create N2 molecule structure
        molecule_structure = molecule('N2')
        molecule_structure.calc = EMT()
        # From ase example:https://wiki.fysik.dtu.dk/ase/tutorials/atomization.html, Obs: they 
        # had 2N Atomization energy(cohesive energy)= 9.76 eV, Which means N alone is 9.76/2 = 4.88 ev/atom
        self.assertTrue((4.82 < calculate_cohesive_energy(molecule_structure)) and
                        (calculate_cohesive_energy(molecule_structure) < 4.93))

    def test_MSD_lindemann_diffusion(self):
        atom1 = FaceCenteredCubic(symbol="Cu", size=(10, 10, 10), pbc=True)

        d = approx_lattice_constant(atom1)
        
        dict1 = {'potential': 'EMT', 'ensemble': 'NVT', 'temperature': 300, 
                 'step_number': 1000, 'time_step': 1, 'friction': 0.005}
        
        value_dict = {'mean_square_displacement': [atom1.get_positions()],
                    'lindemann_criterion': [0],
                    'self_diffusion_coefficient': [0]}
        
        mod_atom = run_simple_md_simulation(atom1, dict1)

        MSD = calc_mean_square_displacement(mod_atom, value_dict)
        lindemann = lindemann_criterion(value_dict, d)
        self_diffusion = self_diffusion_coefficent(value_dict, dict1['step_number']*dict1['time_step'])

        atom2 = FaceCenteredCubic(symbol="Cu", size=(10, 10, 10), pbc=True)
        
        dict2 = {'potential': 'EMT', 'ensemble': 'NVT', 'temperature': 2000, 
                 'step_number': 1000, 'time_step': 1, 'friction': 0.005}
        
        value_dict2 = {'mean_square_displacement': [atom2.get_positions()],
                    'lindemann_criterion': [0],
                    'self_diffusion_coefficient': [0]}
        
        mod_atom2 = run_simple_md_simulation(atom2, dict2)

        MSD2 = calc_mean_square_displacement(mod_atom2, value_dict2)
        lindemann2 = lindemann_criterion(value_dict2, d)
        self_diffusion2 = self_diffusion_coefficent(value_dict2, dict2['step_number']*dict2['time_step'])

        self.assertTrue((MSD < 0.1) and (lindemann < 0.1) and (self_diffusion < 0.001) and
                        (MSD2 > 0.1) and (lindemann2 > 0.1) and (self_diffusion2 > 0.00001))
        

    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
