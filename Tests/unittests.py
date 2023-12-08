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
from Simulation.simple_simulation import run_simple_md_simulation
from Simulation.lattice_constant import optimize_scaling
from Simulation.calc_properties import approx_lattice_constant, calc_temp, calc_pressure, calc_mean_square_displacement, lindemann_criterion, self_diffusion_coefficent
from Simulation.calc_bulk_properties import calc_bulk_modulus, calculate_cohesive_energy, calc_elastic
from Simulation.run_md_simulation import run_single_md_simulation
from Simulation.simple_simulation import run_simple_md_simulation
from Gather_data.download_data import get_ASE_atoms_from_material_id
from User_interface.user_interface import initiate_gui
import numpy as np
from ase import units
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

        # Melting temperature for Copper is 1357.7 K. When using a lower temperature we expect Lindemann < 0.1
        # and low MSD and diffusion. For temperatures greater than 1357.7 we expect Lindemann > 0.1 with greater
        # MSD and diffusion values
        self.assertTrue((MSD < 0.1) and (lindemann < 0.1) and (self_diffusion < 0.001) and
                        (MSD2 > 0.1) and (lindemann2 > 0.1) and (self_diffusion2 > 0.00001))
        

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
        #path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files'
        #traj = Trajectory(path_to_traj_folder + "/NVT_Cu_validation_test.traj", "w")
        #traj.write(atoms)
        
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

        ## Temperature
        # Remove the first 40% of data points since they usually havent reached equilibrium
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
        # Physics handbook give 138 GPa at 300 Kelvin
        avg_bulk = np.average(material_data_dict['bulk_modulus'])
        self.assertTrue(100 < avg_bulk < 170)

        ## Shear modulus
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?chemsys=Cu
        # Physics handbook give 48.3 GPa at 300 Kelvin
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(35 < avg_shear_modulus < 67)
        

        ## Youngs modulus
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        self.assertTrue(100 < avg_youngs_modulus < 200)
        # We get a value of 163 GPa, as compared to 110 on Wiki, but should be alright
        # since its of the same magnitude (discussed with Abijith)

        ## Poisson ratio
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?formula=Cu
        # Expected 0.35
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.1 < avg_poisson_ratio < 0.5)

        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.5)

        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann < 0.2)

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

        ## Temperature
        # Remove the first 40% of data points since they usually havent reached equilibrium
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
        # Taken from Physics Handbook T-1.1, expected 103 GPa
        avg_bulk = np.average(material_data_dict['bulk_modulus'])
        self.assertTrue(63 < avg_bulk < 143)

        ## Shear modulus
        # Taken from Physics Handbook T-1.1, expected 30.3 GPa
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(-10.3 < avg_shear_modulus < 70.3)

        ## Youngs modulus
        # Taken from Physics Handbook T-1.1, expected 82.7 GPa
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        self.assertTrue(42.7 < avg_youngs_modulus < 122.7)

        ## Poisson ratio
        # Taken from https://next-gen.materialsproject.org/materials/mp-124?formula=Ag
        # Expected 0.43
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.2 < avg_poisson_ratio < 0.6)

        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.5)

        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann < 0.2)

        ## Self-diffusion coefficient
        avg_self_diffusion = np.average(material_data_dict['self_diffusion_coefficient'][21:])
        self.assertTrue(avg_self_diffusion < 0.001)

        ## Lattice constant
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/ValidationTests/NVT_Ag_validation_test.traj", 'r')
        atoms = traj[-1]
        nearest_neighbour_distance = approx_lattice_constant(atoms)
        self.assertTrue(2.7 < nearest_neighbour_distance < 3)
        # We get a value of 2.78, that is a lattice constant of 3.94 Å (Expected 4.086 Å)


    def test_simulations_for_validating_code_NVE_Cu(self):
        # In this test requirements 32-35 will all be validated by running 4 large simulations
        # where all the properties are measured. Only the results of the simulations will be 
        # examined here, as the testing would take far too much time otherwise.

        # The following is the atom initialization
        atoms = FaceCenteredCubic(symbol="Cu", size=(10, 10, 10), pbc=True)

        # Location of the trajectory file: Input_trajectory_files/Cu_validation_test.traj
        # Location of configuration file: Input_config_files/NVT_validation test
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/NVE_Cu_validation_test.traj", "w")
        traj.write(atoms)

        #run_single_md_simulation("NVE_validation_test.ini",
        #                         "NVE_Cu_validation_test.traj",
        #                         "ValidationTests/NVE_Cu_validation_test")
        # The simulation results are saved as Cu_validation_test in the ValidationTests
        # folder in Output_traj_files and Output_text_files

        # First all the data is saved in a dictionary
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/ValidationTests/NVE_Cu_validation_test.txt'
        opened_file = open(path, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)

        ## Total energy
        # Remove the first 40% of data points since they usually havent reached equilibrium
        energy_equilibrium = material_data_dict['total_energy'][41:]
        avg_interval = []
        # Check the average temperature in intervals
        for x in range(16):
            avg_interval.append(np.average(energy_equilibrium[x*10:x*10+9]))
        total_avg_temp = np.average(energy_equilibrium)

        # Compare the interval averages to the total average
        diff_temp = total_avg_temp - np.array(avg_interval)
        #self.assertTrue(abs(np.average(diff_temp)) < 10)   

        ## Pressure
        avg_pressure = np.average(material_data_dict['pressure'][41:])
        self.assertTrue(abs(avg_pressure) < 0.5)

        ## Bulk modulus
        # Taken from Physics Handbook T-1.1, expected 151 GPa
        avg_bulk = np.average(material_data_dict['bulk_modulus'])
        self.assertTrue(121 < avg_bulk < 181)

        ## Shear modulus
        # Taken from Physics Handbook T-1.1, expected 57 GPa
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(27 < avg_shear_modulus < 87)

        ## Youngs modulus
        # Taken from Physics Handbook T-1.1, expected 110 GPa
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        self.assertTrue(60 < avg_youngs_modulus < 180)

        ## Poisson ratio
        # Taken from https://next-gen.materialsproject.org/materials/mp-30?formula=Cu
        # Expected 0.35
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.1 < avg_poisson_ratio < 0.5)

        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.2)

        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann < 0.2)

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
        # It makes sense that we underestimate it slightly if the atom has more than 4 nearest neighboors
        # (since we will pick the ones that happens to be closest at the specific moment)


    def test_simulations_for_validating_code_NVE_Au(self):
        # In this test requirements 32-35 will all be validated by running 4 large simulations
        # where all the properties are measured. Only the results of the simulations will be 
        # examined here, as the testing would take far too much time otherwise.

        # The following is the atom initialization
        atoms = FaceCenteredCubic(symbol="Au", size=(10, 10, 10), pbc=True)

        # Location of the trajectory file: Input_trajectory_files/Cu_validation_test.traj
        # Location of configuration file: Input_config_files/NVT_validation test
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/NVE_Au_validation_test.traj", "w")
        traj.write(atoms)

        #run_single_md_simulation("NVE_validation_test.ini",
        #                         "NVE_Au_validation_test.traj",
        #                         "ValidationTests/NVE_Au_validation_test")
        # The simulation results are saved as Cu_validation_test in the ValidationTests
        # folder in Output_traj_files and Output_text_files

        # First all the data is saved in a dictionary
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/ValidationTests/NVE_Au_validation_test.txt'
        opened_file = open(path, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)

        ## Total energy
        # Remove the first 40% of data points since they usually havent reached equilibrium
        energy_equilibrium = material_data_dict['total_energy'][41:]
        avg_interval = []
        # Check the average temperature in intervals
        for x in range(16):
            avg_interval.append(np.average(energy_equilibrium[x*10:x*10+9]))
        total_avg_energy = np.average(energy_equilibrium)

        # Compare the interval averages to the total average
        diff_energy = total_avg_energy - np.array(avg_interval)
        self.assertTrue(abs(np.average(diff_energy)) < 10)

        ## Pressure
        avg_pressure = np.average(material_data_dict['pressure'][41:])
        self.assertTrue(abs(avg_pressure) < 0.5)

        ## Bulk modulus
        # Taken from Physics Handbook T-1.1, expected 217 GPa
        avg_bulk = np.average(material_data_dict['bulk_modulus'])
        self.assertTrue(137 < avg_bulk < 297)

        ## Shear modulus
        # Taken from Physics Handbook T-1.1, expected 27 GPa
        avg_shear_modulus = np.average(material_data_dict['shear_modulus'])
        self.assertTrue(-17 < avg_shear_modulus < 67)

        ## Youngs modulus
        # Taken from Physics Handbook T-1.1, expected 78.5 GPa
        avg_youngs_modulus = np.average(material_data_dict['youngs_modulus'])
        self.assertTrue(48.5 < avg_youngs_modulus < 118.5)

        ## Poisson ratio
        # Taken from https://next-gen.materialsproject.org/materials/mp-81?formula=Au
        # Expected 0.52
        avg_poisson_ratio = np.average(material_data_dict['poisson_ratio'])
        self.assertTrue(0.32 < avg_poisson_ratio < 0.72)

        ## MSD
        avg_MSD = np.average(material_data_dict['mean_square_displacement'][21:])
        self.assertTrue(avg_MSD < 0.3)
        # This is expected to be constantly low when the simulation temperature is far below 
        # the materials melting point

        ## Lindemann criterion
        avg_lindemann = np.average(material_data_dict['lindemann_criterion'][21:])
        self.assertTrue(avg_lindemann < 0.2)

        ## Self-diffusion coefficient
        avg_self_diffusion = np.average(material_data_dict['self_diffusion_coefficient'][21:])
        self.assertTrue(avg_self_diffusion < 0.001)

        ## Lattice constant
        path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files'
        traj = Trajectory(path_to_traj_folder + "/ValidationTests/NVE_Au_validation_test.traj", 'r')
        atoms = traj[-1]
        nearest_neighbour_distance = approx_lattice_constant(atoms)
        self.assertTrue(2.5 < nearest_neighbour_distance < 2.9)
        # We get a value of 2.77, that is a lattice constant of 3.92 Å (Expected 4.078 Å)


    def test_GUI(self):
        # There will be further testing when other methods connected to the gui has been developed.
        gui = initiate_gui()
        self.assertTrue(type(gui) == Tk)

    def test_specific_heat_capacity(self):
        # Get the directory path
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Premade_simulation_data/Specific_heat_capacity/'

        # Get a list of all files in the directory
        all_files = os.listdir(path)

        for file_name in all_files:
            # NVE ensemble tests
            # Copper
            if "NVE_copper_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    #print("spec heat capa: ", specific_heat_capacity)

                    # Check if the value is within the specified range
                    # The experimental value from Physics handbook T-1.1 at 300K for copper = 385 J/Kg*K
                    self.assertTrue(383 <= specific_heat_capacity <= 393)

            # Silver
            if "NVE_silver_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    # print("spec heat capa: ", specific_heat_capacity)

                    # The experimental value from Physics handbook T-1.1 at 300K for silver = 235 J/Kg*K
                    self.assertTrue(230 <= specific_heat_capacity <= 240)

            # Aluminium
            if "NVE_aluminium_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    # print("spec heat capa: ", specific_heat_capacity)

                    # The experimental value from Physics handbook T-1.1 at 300K for Aluminium = 897 J/Kg*K
                    self.assertTrue(895 <= specific_heat_capacity <= 930)

            # NVT ensemble tests
            # Copper (Better value than NVE ensemble)
            if "NVT_copper_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    # print("spec heat capa: ", specific_heat_capacity)

                    # Check if the value is within the specified range
                    # The experimental value from Physics handbook T-1.1 at 300K for copper = 385 J/Kg*K
                    self.assertTrue(383 <= specific_heat_capacity <= 393)
            
            # Silver (worse value than NVE ensemble)
            if "NVT_silver_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    # print("spec heat capa: ", specific_heat_capacity)

                    # The experimental value from Physics handbook T-1.1 at 300K for silver = 235 J/Kg*K
                    self.assertTrue(200 <= specific_heat_capacity <= 240)

            # Aluminium (So much worse value than NVE ensemble))
            if "NVT_aluminium_heat_capa_3timestep.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract specific_heat_capacity from the data
                    specific_heat_capacity = data.get("specific_heat_capacity", [])[0]
                    # print("spec heat capa: ", specific_heat_capacity)

                    # The experimental value from Physics handbook T-1.1 at 300K for Aluminium = 897 J/Kg*K
                    self.assertTrue(730 <= specific_heat_capacity <= 900)

    def test_debye_temperature(self):
        # Get the directory path
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Premade_simulation_data/Debye_temerpature/'

        # Get a list of all files in the directory
        all_files = os.listdir(path)

        for file_name in all_files:
            # NVE ensemble tests
            # Copper
            if "NVE_copper_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for copper = 343 K
                    self.assertTrue(330 <= debye_temperature <= 500)

            # Silver
            if "NVE_silver_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for silver = 225 K
                    self.assertTrue(220 <= debye_temperature <= 350)

            # Aluminium
            if "NVE_aluminium_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for Aluminium = 428 K
                    self.assertTrue(420 <= debye_temperature <= 450)

            # NVT ensemble tests
            # Copper
            if "NVT_copper_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for copper = 343 K
                    self.assertTrue(330 <= debye_temperature <= 500)

            # Silver
            if "NVT_silver_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for silver = 225 K
                    self.assertTrue(220 <= debye_temperature <= 350)

            # Aluminium
            if "NVT_aluminium_debye_temp.txt" in file_name:
                file_path = os.path.join(path, file_name)
                with open(file_path, 'r') as txt_file:
                    # Load JSON data from the file
                    data = json.load(txt_file)

                    # Extract debye temperature from the data
                    debye_temperature = data.get("time_average_of_debye_temperature", [])
                    #print("debye temp: ", debye_temperature)

                    # Check if the value is within the specified range
                    # The experimental value from "Introduction to Solid State Physics" by Charles Kittel page 116 at 300K for Aluminium = 428 K
                    self.assertTrue(420 <= debye_temperature <= 450)


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
    
    
    def test_lattice_constant(self):
        # Lattice constant for Cu (fcc) is 3.61 Å, which gives a nearest neighbor distance of 2.55 Å
        atom_Cu = FaceCenteredCubic(symbol="Cu", size=(5, 5, 5), pbc=True)
        self.assertTrue(2.4 <= approx_lattice_constant(atom_Cu) <= 2.7)
        # Lattice constant for Ag (fcc) is 4.09 Å, which gives a nearest neighbor distance of 2.89 Å
        atom_Ag = FaceCenteredCubic(symbol="Ag", size=(5, 5, 5), pbc=True)
        self.assertTrue(2.8 <= approx_lattice_constant(atom_Ag) <= 3.05)
        
        dict1 = {'potential': 'EMT', 'ensemble': 'NVT', 'temperature': 300, 
                 'step_number': 2000, 'time_step': 1, 'friction': 0.005}
        
        value_dict1 = {'mean_square_displacement': [atom_Cu.get_positions()],
                    'lindemann_criterion': [0],
                    'self_diffusion_coefficient': [0]}
        
        mod_atom_Cu = run_simple_md_simulation(atom_Cu, dict1)
        self.assertTrue(2.3 <= approx_lattice_constant(mod_atom_Cu) <= 2.7)
        mod_atom_Ag = run_simple_md_simulation(atom_Ag, dict1)
        self.assertTrue(2.7 <= approx_lattice_constant(mod_atom_Ag) <= 3.05)
        
        # The nearest neighbour distance varied between 2.48-2.39 for Copper and
        # between 2.81-2.70 for Gold when varying temperature between 300-1900 K.
        # That is lattice constant for Cu: 3.50-3.37 Å and Ag: 3.97-3.81 Å.



if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(UnitTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
