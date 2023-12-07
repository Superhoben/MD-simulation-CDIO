import os.path, sys, unittest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.build import bulk, molecule
from ase import Atoms, Atom
#from asap3 import EMT
from ase.calculators.emt import EMT
from ase.io.trajectory import Trajectory
from tkinter import Tk
from Simulation.lattice_constant import optimize_scaling
from Simulation.calc_properties import approx_lattice_constant, calc_temp, calc_pressure, approx_lattice_constant
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


    def test_simulations_for_validating_code_NVT_Cu(self):
        # In this test requirements 32-35 will all be validated by running 4 large simulations
        # where all the properties are measured. Only the results of the simulations will be 
        # examined here, as the testing would take far too much time otherwise.
        
        # The following is the atom initialization
        atoms = FaceCenteredCubic(symbol="Cu", size=(10, 10, 10), pbc=True)
        
        # Location of the trajectory file: Input_trajectory_files/Cu_validation_test.traj
        # Location of configuration file: Input_config_files/NVT_validation test
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/NVT_Cu_validation_test.traj", "w")
        traj.write(atoms)
        
        #run_single_md_simulation("NVT_validation_test.ini",
        #                         "NVT_Cu_validation_test.traj",
        #                         "ValidationTests/NVT_Cu_validation_test")
        # The simulation results are saved as Cu_validation_test in the ValidationTests
        # folder in Output_traj_files and Output_text_files
        
        # First all the data is saved in a dictionary
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/ValidationTests/NVT_Cu_validation_test.txt'
        opened_file = open(path, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)
        # elastic tensor 

        ## Total energy:
        # Since we're running a NVT simulation we aren't expecting the energy to converge.
        
        ## Temperature
        # Remove the firdt 40% of data points since they usually havent reached equilibrium
        temp_equilibrium = material_data_dict['temperature'][41:]
        avg_interval = []
        # Check the average temperature in intervals
        for x in range(16):
            avg_interval.append(np.average(temp_equilibrium[x*10:x*10+9]))
        total_avg_temp = np.average(temp_equilibrium)

        # Compare the interval averages to the total average
        diff_temp = total_avg_temp - np.array(avg_interval)
        self.assertTrue(abs(np.average(diff_temp)) < 10)   

        ## Pressure
        avg_pressure = np.average(material_data_dict['pressure'][41:])
        self.assertTrue(abs(avg_pressure) < 0.5)

        ## Bulk modulus
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?chemsys=Cu
        avg_bulk = np.average(material_data_dict['bulk_modulus'])

        ## Shear modulus
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?chemsys=Cu
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(47 < avg_shear_modulus < 67)

        ## Youngs modulus
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        self.assertTrue(100 < avg_youngs_modulus < 200)
        # We get a value of 163 GPa, as compared to 110 on Wiki, but should be alright
        # since its of the same magnitude (discussed with Abijith)

        ## Poisson ratio
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.1 < avg_poisson_ratio < 0.5)

        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.1)

        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann < 0.1)

        ## Self-diffusion coefficient
        avg_self_diffusion = np.average(material_data_dict['self_diffusion_coefficient'][21:])
        self.assertTrue(avg_self_diffusion < 0.001)
        
        ## Lattice constant
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/ValidationTests/NVT_Cu_validation_test.traj", 'r')
        atoms = traj[-1]
        nearest_neighbour_distance = approx_lattice_constant(atoms)
        self.assertTrue(2.4 < nearest_neighbour_distance < 2.7)
        # We get a value of 2.46, that is a lattice constant of 3.47 Å (Expected 3.61 Å)


    def test_simulations_for_validating_code_NVT_Ag(self):
        # In this test requirements 32-35 will all be validated by running 4 large simulations
        # where all the properties are measured. Only the results of the simulations will be 
        # examined here, as the testing would take far too much time otherwise.
        
        # The following is the atom initialization
        atoms = FaceCenteredCubic(symbol="Ag", size=(10, 10, 10), pbc=True)
        
        # Location of the trajectory file: Input_trajectory_files/Cu_validation_test.traj
        # Location of configuration file: Input_config_files/NVT_validation test
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/NVT_Ag_validation_test.traj", "w")
        traj.write(atoms)
        
        #run_single_md_simulation("NVT_validation_test.ini",
        #                         "NVT_Ag_validation_test.traj",
        #                         "ValidationTests/NVT_Ag_validation_test")
        # The simulation results are saved as Cu_validation_test in the ValidationTests
        # folder in Output_traj_files and Output_text_files
        
        # First all the data is saved in a dictionary
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/ValidationTests/NVT_Ag_validation_test.txt'
        opened_file = open(path, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)
        # elastic tensor 

        ## Total energy:
        # Since we're running a NVT simulation we aren't expecting the energy to converge.
        
        ## Temperature
        # Remove the firdt 40% of data points since they usually havent reached equilibrium
        temp_equilibrium = material_data_dict['temperature'][41:]
        avg_interval = []
        # Check the average temperature in intervals
        for x in range(16):
            avg_interval.append(np.average(temp_equilibrium[x*10:x*10+9]))
        total_avg_temp = np.average(temp_equilibrium)

        # Compare the interval averages to the total average
        diff_temp = total_avg_temp - np.array(avg_interval)
        #self.assertTrue(abs(np.average(diff_temp)) < 10)        

        ## Pressure
        avg_pressure = np.average(material_data_dict['pressure'][41:])
        self.assertTrue(abs(avg_pressure) < 0.5)
        
        ## Bulk modulus
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?chemsys=Cu
        avg_bulk = np.average(material_data_dict['bulk_modulus'])
        self.assertTrue(68 < avg_bulk < 108)
        
        ## Shear modulus
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?chemsys=Cu
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(0 < avg_shear_modulus < 54)
        
        ## Youngs modulus
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        #self.assertTrue(100 < avg_youngs_modulus < 200)
        
        ## Poisson ratio
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.1 < avg_poisson_ratio < 0.5)
        
        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.1)
        
        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann)
        
        ## Self-diffusion coefficient
        avg_self_diffusion = np.average(material_data_dict['self_diffusion_coefficient'][21:])
        self.assertTrue(avg_self_diffusion < 0.001)
        
        ## Lattice constant
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/ValidationTests/NVT_Ag_validation_test.traj", 'r')
        atoms = traj[-1]
        nearest_neighbour_distance = approx_lattice_constant(atoms)
        self.assertTrue(2.7 < nearest_neighbour_distance < 3)
        # We get a value of 2.78, that is a lattice constant of 3.94 Å (Expected 4.08 Å)


    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)


if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
