import os.path, sys, unittest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.build import bulk, molecule
from ase import Atoms, Atom
#from asap3 import EMT
from ase.calculators.emt import EMT
from tkinter import Tk
from Simulation.lattice_constant import optimize_scaling
from Simulation.calc_properties import calc_temp, calc_pressure, approx_lattice_constant
from Simulation.calc_bulk_properties import calc_bulk_modulus, calculate_cohesive_energy, calc_elastic
from Simulation.run_md_simulation import run_single_md_simulation
from Gather_data.download_data import get_ASE_atoms_from_material_id
from User_interface.user_interface import initiate_gui
from API_key import start_program


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

    def test_lattice_constant_gradient_descent(self):
        """Uses a material with known lattice constant, disorts it and see if the lattice method is able
        to find the original lattice constant which also should be the optimal one."""
        #atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(2, 2, 2),
                          pbc=True)
        atoms.set_cell(atoms.cell*0.5, scale_atoms=True)  # Rescale unit cell so atoms are now in suboptimal positions
        atoms.calc = EMT()
        scaling = optimize_scaling(atoms, {'optimal_scaling': [], 'iterations_to_find_scaling': []})
        # print(scaling)
        self.assertTrue((0.98 < scaling*0.5) and (scaling*0.5 < 1.02))
        

    # More test are needed for this function
    def test_calc_pressure_no_field(self):
        #atoms = get_ASE_atoms_from_material_id('mp-30')  # Get atoms object with atoms in optimal positions
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(2, 2, 2),
                          pbc=True)
        atoms.calc = EMT()  # If atoms are only, Ni, Cu, Pd, Ag, Pt or Au no input parameters are nessecary
        # When atoms are still at optimal positions no pressure should occur, +-0.1 GPa tolerance is given
        self.assertTrue((-0.1 < calc_pressure(atoms)) and (calc_pressure(atoms) < 0.1))
        
    def test_bulk_modulus(self):
        # lattice_constant = 4
        # choosing more specific silver to test the bulk modulus
        #atoms = get_ASE_atoms_from_material_id('mp-124')
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Ag",
                          size=(2, 2, 2),
                          pbc=True)
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
                        (calculate_cohesive_energy(atom_structure_Ag) < 3))
        self.assertTrue((3.76 < calculate_cohesive_energy(atom_structure_Au)) and
                        (calculate_cohesive_energy(atom_structure_Au) < 3.86))
        self.assertTrue((4.39 < calculate_cohesive_energy(atom_structure_Ni)) and
                        (calculate_cohesive_energy(atom_structure_Ni) < 4.49))
        self.assertTrue((3.44 < calculate_cohesive_energy(atom_structure_Cu)) and
                        (calculate_cohesive_energy(atom_structure_Cu) < 3.54))

    def test_cohesive_energy(self):
        # Create N2 molecule structure
        molecule_structure = molecule('N2')
        molecule_structure.calc = EMT()
        # From ase example:https://wiki.fysik.dtu.dk/ase/tutorials/atomization.html, Obs: they 
        # had 2N Atomization energy(cohesive energy)= 9.76 eV, Which means N alone is 9.76/2 = 4.88 ev/atom
        self.assertTrue((4.82 < calculate_cohesive_energy(molecule_structure)) and
                        (calculate_cohesive_energy(molecule_structure) < 4.93))

    def test_approx_lattice(self):
        # Lattice constant for Cu (fcc) is 3.61 Å, which gives a nearest neighbor distance of 2.55 Å
        atoms = FaceCenteredCubic(symbol="Cu", size=(5, 5, 5), pbc=True)
        self.assertTrue(2.4 <= approx_lattice_constant(atoms) <= 2.7)
        # Lattice constant for Ag (fcc) is 4.09 Å, which gives a nearest neighbor distance of 2.89 Å
        atoms = FaceCenteredCubic(symbol="Ag", size=(5, 5, 5), pbc=True)
        self.assertTrue(2.8 <= approx_lattice_constant(atoms) <= 3.05)

    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)

    def test_elastic(self):
        # Elastic properties for different materials.
        # Young's- Bulk- and Shear modulus from Physics Handbook T-1.1
        # Poisson's ratio from https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html (can anyone find a better source?)
        # Element | Young's modulus (GPa) | Bulk modulus (GPa) | Shear modulus (GPa) | Poisson's ratio |
        # Cu      |        129.8          |        137.8       |        48.3         |       0.36      |
        # Ag      |         82.7          |        103.6       |        30.3         |       0.37      |
        # Pt      |          168          |          228       |          61         |        0.3      |

        # Cu
        atoms_Cu = FaceCenteredCubic(symbol="Cu", size=(2, 2, 2), pbc=True, latticeconstant = 3.61)
        atoms_Cu.calc = EMT()
        Cij_Cu, Bij_Cu, cij_order_Cu, bulk_modulus_Cu, shear_modulus_Cu, \
        youngs_modulus_Cu, poisson_ratio_Cu = calc_elastic(atoms_Cu)

        self.assertTrue(129.8*0.8 < youngs_modulus_Cu < 129.8*1.2)
        self.assertTrue(137.8*0.8 < bulk_modulus_Cu < 137.8*1.2)
        self.assertTrue(48.3*0.75 < shear_modulus_Cu < 48.3*1.25)
        self.assertTrue(0.36*0.8 < poisson_ratio_Cu < 0.36*1.2)

        # Ag
        atoms_Ag = FaceCenteredCubic(symbol="Ag", size=(2, 2, 2), pbc=True)
        atoms_Ag.calc = EMT()
        Cij_Ag, Bij_Ag, cij_order_Ag, bulk_modulus_Ag, shear_modulus_Ag, \
        youngs_modulus_Ag, poisson_ratio_Ag = calc_elastic(atoms_Ag)

        self.assertTrue(82.7*0.8 < youngs_modulus_Ag < 82.7*1.2)
        self.assertTrue(103.6*0.8 < bulk_modulus_Ag < 103.6*1.2)
        self.assertTrue(30.3*0.75 < shear_modulus_Ag < 30.3*1.25)
        self.assertTrue(0.37*0.9 < poisson_ratio_Ag < 0.37*1.1)

        # Pt
        atoms_Pt = FaceCenteredCubic(symbol="Pt", size=(2, 2, 2), pbc=True)
        atoms_Pt.calc = EMT()
        Cij_Pt, Bij_Pt, cij_order_Pt, bulk_modulus_Pt, shear_modulus_Pt, \
        youngs_modulus_Pt, poisson_ratio_Pt = calc_elastic(atoms_Pt)

        self.assertTrue(168*0.8 < youngs_modulus_Pt < 168*1.2)
        self.assertTrue(228*0.5 < bulk_modulus_Pt < 228*1.5)
        self.assertTrue(61*0.5 < shear_modulus_Pt < 61*1.5)
        self.assertTrue(0.3*0.6 < poisson_ratio_Pt < 0.3*1.4)


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
