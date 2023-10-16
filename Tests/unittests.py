import os, sys, unittest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from Simulation import *
from Gather_data import *
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3 import EMT


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
        self.assertTrue(270 < calc_properties.calc_temp(atoms) < 330)

    def test_lattice_constant(self):
        lattice_const = calc_properties.optimize_lattice_const(
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
        atoms = gather_data.get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms.set_cell(atoms.cell*0.5, scale_atoms=True)  # Rescale unit cell so atoms are now in suboptimal positions
        scaling, _, _ = lattice_constant.optimize_lattice_const_gradient_descent(atoms, lattice_constant.example_simulation_function)
        print(scaling)
        self.assertTrue((0.98 < scaling*0.5) and (scaling*0.5 < 1.02))

    # More test are needed for this function
    def test_calc_pressure_no_field(self):
        atoms = gather_data.get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms.calc = EMT()  # If atoms are only, Ni, Cu, Pd, Ag, Pt or Au no input parameters are nessecary
        # When atoms are still at optimal positions no pressure should occur, +-0.1 GPa tolerance is given
        self.assertTrue((-0.1 < calc_properties.calc_pressure(atoms)) and (calc_properties.calc_pressure(atoms) < 0.1))

    def test_bulk_modulus(self):
        # lattice_constant = 4
        # choosing more specific silver to test the bulk modulus
        atoms = get_ASE_atoms_from_material_id('mp-124')
        atoms.calc = EMT()
        create_traj_file(atoms,lattice_constant = 4)
        # Read the traj file from the first atom to the 10th atom
        calc_bulk_modulus("atoms.traj@0:9")
        # From material website https://next-gen.materialsproject.org/materials/mp-124?chemsys=Ag#elastic_constants
        # we have the bulk modulus for Ag to be 88 GPa
        # According to Rickard, it is absolutely no problem to get 100 GPa of bulk modulus since we are first taking the
        # second derivative approximation and secondly we are using EMT calculator.

        self.assertTrue((86 < calc_bulk_properties.calc_bulk_modulus("atoms.traj@0:9")) and (calc_bulk_properties.calc_bulk_modulus("atoms.traj@0:9") < 105))


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
