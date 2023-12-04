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
from Simulation.calc_bulk_properties import calc_bulk_modulus, calculate_cohesive_energy
from Simulation.run_md_simulation import run_single_md_simulation
from Gather_data.download_data import get_ASE_atoms_from_material_id
from User_interface.user_interface import initiate_gui
from API_key import start_program
import json
import numpy as np


opened_file_NVE = open("Tests/SimulationOutputs/NVE_300K_for_testing.txt", "r")
opened_file_NVT = open("Tests/SimulationOutputs/NVT_300K_for_testing.txt", "r")
data_NVE = opened_file_NVE.readline()
data_NVT = opened_file_NVT.readline()
opened_file_NVE.close()
opened_file_NVT.close()
material_data_dict_NVE = json.loads(data_NVE)
material_data_dict_NVT = json.loads(data_NVT)


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


    def test_cohesive_energy2(self):
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


    def test_internal_pressure(self):
        # Remove the first points since theyre usually havent reached equilibrium
        # For an optimized volume the internal presssure should be minimized
        avg_pressure_NVE = np.average(material_data_dict_NVE['pressure'][11:])
        avg_pressure_NVT = np.average(material_data_dict_NVT['pressure'][11:])
        
        self.assertTrue(abs(avg_pressure_NVE) < 0.5)
        self.assertTrue(abs(avg_pressure_NVT) < 0.5)


    def test_reached_equilibrium(self):
        # Remove the first points since theyre usually havent reached equilibrium
        NVE_equilibrium = material_data_dict_NVE['total_energy'][11:]
        NVT_equilibrium = material_data_dict_NVT['temperature'][11:]
        avg_interval_NVE = []
        avg_interval_NVT = []
        # Check the average energy/temperature in intervals
        for x in range(9):
            avg_interval_NVE.append(np.average(NVE_equilibrium[x*10:x*10+9]))
            avg_interval_NVT.append(np.average(NVT_equilibrium[x*10:x*10+9]))
        total_avg_NVE = np.average(NVE_equilibrium)
        total_avg_NVT = np.average(NVT_equilibrium)
        
        # Compare the interval averages to the total average
        diff_NVE = total_avg_NVE - np.array(avg_interval_NVE)
        diff_NVT = total_avg_NVT - np.array(avg_interval_NVT)

        self.assertTrue(abs(np.average(diff_NVE)) < 0.001)
        # Since the temperature differs more than the energy for NVT 
        # a bigger error margin is used
        self.assertTrue(abs(np.average(diff_NVT)) < 10)


    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
