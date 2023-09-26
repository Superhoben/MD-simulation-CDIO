import sys, unittest
from firstunittest import *

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
        

if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
