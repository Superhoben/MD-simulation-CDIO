import sys, unittest
from firstunittest import *
from calc_properties import calc_temp
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


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
        self.assertTrue(290 < calc_temp(atoms) < 310)


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
