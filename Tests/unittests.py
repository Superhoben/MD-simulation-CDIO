import os.path, sys, unittest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.build import bulk, molecule
from ase import Atoms
#from asap3 import EMT
from ase.calculators.emt import EMT
from tkinter import Tk
from Simulation.lattice_constant import optimize_lattice_const, example_simulation_function, optimize_lattice_const_gradient_descent
from Simulation.calc_properties import calc_temp, calc_pressure
from Simulation.calc_bulk_properties import create_traj_file, calc_bulk_modulus, calculate_cohesive_energy
from Simulation.run_md_simulation import run_NVE_NVT
from Gather_data.download_data import get_ASE_atoms_from_material_id
from User_interface.user_interface import initiate_gui


class UnitTests(unittest.TestCase):
    """Definitions of unittests.

    Different tests can be added here from different files. To add a new test,
    import the functions from the file that is to be tested and then write a 
    new test function below the existing ones.

    See example tests below.
    """
    def test_calc_temp(self):
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(5, 5, 5),
                          pbc=True)
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)
        self.assertTrue(270 < calc_temp(atoms) < 330)

    def test_lattice_constant(self):
        lattice_const = optimize_lattice_const(
           FaceCenteredCubic(
               directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]],
               size=(2, 2, 3),
               symbol="Cu",
               pbc=(1, 1, 0),
                )
            )
        self.assertTrue((0.95 < lattice_const) and (lattice_const < 1.05))

    def test_lattice_constant_gradient_descent(self):
        """Uses a material with known lattice constant, disorts it and see if the lattice method is able
        to find the original lattice constant which also should be the optimal one."""
        atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms.set_cell(atoms.cell*0.5, scale_atoms=True)  # Rescale unit cell so atoms are now in suboptimal positions
        scaling, _, _ = optimize_lattice_const_gradient_descent(atoms, example_simulation_function)
        self.assertTrue((0.98 < scaling*0.5) and (scaling*0.5 < 1.02))

    # More test are needed for this function
    def test_calc_pressure_no_field(self):
        atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms.calc = EMT()  # If atoms are only, Ni, Cu, Pd, Ag, Pt or Au no input parameters are nessecary
        # When atoms are still at optimal positions no pressure should occur, +-0.1 GPa tolerance is given
        self.assertTrue((-0.1 < calc_pressure(atoms)) and (calc_pressure(atoms) < 0.1))

    def test_bulk_modulus(self):
        # lattice_constant = 4
        # choosing more specific silver to test the bulk modulus
        atoms = get_ASE_atoms_from_material_id('mp-124')
        atoms.calc = EMT()
        # Another way to create our object "bulk" this time
        #ag_bulk = bulk(name= "Ag",crystalstructure= "fcc", a=4.09)
        #ag_bulk.calc = EMT
        create_traj_file(atoms, lattice_constant=4)
        # Read the traj file from the first atom to the 10th atom
        calc_bulk_modulus("atoms.traj@0:9")
        # From material website https://next-gen.materialsproject.org/materials/mp-124?chemsys=Ag#elastic_constants
        # we have the bulk modulus for Ag to be 88 GPa
        # According to Rickard, it is absolutely no problem to get 100 GPa of bulk modulus since we are first taking the
        # second derivative approximation and secondly we are using EMT calculator.
        self.assertTrue((86 < calc_bulk_modulus("atoms.traj@0:9")) and (calc_bulk_modulus("atoms.traj@0:9") < 105))

    def test_cohesive_energy(self):
        # Create an isolated Ag atom object
        atom_structure = Atoms("Ag", positions=[(0, 0, 0)])
        atom_structure.calc = EMT()
        # Create Ag bulk crystal structure
        bulk_structure = bulk(name="Ag", crystalstructure="fcc", a = 4.09)
        bulk_structure.calc = EMT()
        # From Charles Kittle book "Introduction to Solid State Physics" page 50
        # you can find every single cohesive energy per atom for each elements in eV/atom
        #expected_cohesive_energy_range = (2.8, 3)
        #cohesive_energy = calculate_cohesive_energy(atom_structure, bulk_structure)
        #self.assertTrue(expected_cohesive_energy_range[0] <= cohesive_energy <= expected_cohesive_energy_range[1])
        # Another way to test stuff (1 line of code)
        self.assertTrue((2.8 < calculate_cohesive_energy(atom_structure, bulk_structure)) and
                        (calculate_cohesive_energy(atom_structure,bulk_structure) < 3))

    def test_cohesive_energy(self):
        atom_structure = Atoms("N")
        atom_structure.calc = EMT()
        # Create 2N molecule structure
        #molecule_structure = molecule('N2')
        molecule_structure = Atoms('2N', [(0., 0., 0.), (0., 0., 1.1)])
        molecule_structure.calc = EMT()
        # From ase example:https://wiki.fysik.dtu.dk/ase/tutorials/atomization.html
        self.assertTrue((4.7 < calculate_cohesive_energy(atom_structure, molecule_structure)) and
                        (calculate_cohesive_energy(atom_structure, molecule_structure) < 5))

    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)

    def test_run_NVE_NVT(self):
        atoms = FaceCenteredCubic(size=(7, 7, 7), symbol='Cu', pbc=False)
        config_data = {"show_properties": False, "calc_properties": False, "temperature": 300, "time_step": 5, "interval": 1000, "iterations": 500, "potential": 'EMT', "friction": None}
        new_atoms = run_NVE_NVT(atoms, config_data, 'NVT')
        self.assertTrue(280<calc_temp(new_atoms)<320)

if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
