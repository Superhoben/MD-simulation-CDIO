import sys, unittest
from firstunittest import *
from calc_properties import calc_temp, calc_pressure
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from gather_data import get_ASE_atoms_from_material_id
from asap3 import EMT
from lattice_constant import optimize_lattice_const


class UnitTests(unittest.TestCase):
    """Definitions of unittests.
    
    Different tests can be added here from different files. To add a new test,
    import the functions from the file that is to be tested and then write a 
    new test function below the existing ones.
    
    See example tests below.
    """
    def test1(self):
        check_test = 1
        self.assertTrue(check_test == test1())
        
    def test2(self):
        fib_num_1 = 1
        fib_num_5 = 5
        fib_num_12 = 144
        self.assertTrue(fib_num_1 == calculate_fibonacci_number(1) and
                        fib_num_5 == calculate_fibonacci_number(5) and
                        fib_num_12 == calculate_fibonacci_number(12))
        
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
        self.assertTrue(0.95<lattice_const<1.05)

    # More test are needed for this function
    def test_calc_pressure_no_field(self):
        atoms = get_ASE_atoms_from_material_id('mp-30') # Get atoms object with atoms in optimal positions
        atoms.calc = EMT() # If atoms are only, Ni, Cu, Pd, Ag, Pt or Au no input parameters are nessecary
        # When atoms are still at optimal positions no pressure should occur, +-0.1 GPa tolerance is given
        self.assertTrue(-0.1 < calc_pressure(atoms) < 0.1)
        print(calc_pressure(atoms))


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
