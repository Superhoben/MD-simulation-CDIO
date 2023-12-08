import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from Simulation import simple_simulation
from ase.lattice.cubic import FaceCenteredCubic
from asap3 import EMT


def optimize_scaling(atoms, output_dict={'optimal_scaling': [], 'iterations_to_find_scaling': []}, learning_rate=0.05):
    """Find the optimal scaling of the lattice constant using a gradient descent method.

    The gradient descent method is modified by a sigmoid function since energy changes
    can become very extreme as some potential depend on the -12th power of the distance
    so starting with half the correct distance can result in approximately 4000 times the
    energy which would make a normal gradient search go crazy even with a reasonable
    learning rate. The new scaling s_(k+1) is calculated as:
    e_gradient = (e_(k) - e_(k-1)) / (s_(k) - s_(k-1))
    s_(k+1) = s_(k) - learning_rate*sigmoid(e_gradient)
    The energy at a certain scaling is calculated by using the simulation_function.

    Args:
        atoms(ASE atoms object): The configuration to find an optimal lattice constant for
        output_dict(dict): Dictionary to append the result to.
        learning_rate(float): The scaling converges quicker for larger values but if
             it's too large the scaling will oscillate and not converge at all

    Returns:
        scaling(float): Return the scaling factor which would give the inputed atoms object
            the lowest possible energy
    """
    old_scaling = 1
    scaling = 1.001
    e_scaling_gradient = 0
    number_of_iterations = 0
    old_energy_per_atom = atoms.get_total_energy()
    # Improve the lattice constant performing gradient descent until the gradient becomes
    # sufficiently small. A maximum of 100
    while (abs(e_scaling_gradient) > 0.01) or (number_of_iterations < 3):
        atoms_scaled = atoms.copy()
        atoms_scaled.calc = atoms.calc
        atoms_scaled.set_cell(atoms.cell*scaling, scale_atoms=True)
        energy_per_atom = atoms_scaled.get_total_energy()/len(atoms_scaled)
        e_scaling_gradient = (energy_per_atom-old_energy_per_atom)/(scaling-old_scaling)
        old_scaling = scaling
        scaling = old_scaling - learning_rate*e_scaling_gradient/(2+abs(e_scaling_gradient))
#        print("Scaling gradient: ", e_scaling_gradient)
#        print("Energy per atom: ", energy_per_atom)
#        print("New scaling: ", scaling)
        old_energy_per_atom = energy_per_atom
        number_of_iterations += 1
        if (number_of_iterations > 100):
            scaling = 0
            break
    output_dict['optimal_scaling'].append(scaling)
    output_dict['iterations_to_find_scaling'].append(number_of_iterations)
    return scaling


def optimize_scaling_using_simulation(atoms, simulation_settings: dict,
                                      output_dict={'optimal_scaling': [], 'iterations_to_find_scaling': []}, learning_rate=0.05):
    """Find the optimal scaling of the lattice constant using a gradient descent method.

    The gradient descent method is modified by a sigmoid function since energy changes
    can become very extreme as some potential depend on the -12th power of the distance
    so starting with half the correct distance can result in approximately 4000 times the
    energy which would make a normal gradient search go crazy even with a reasonable
    learning rate. The new scaling s_(k+1) is calculated as:
    e_gradient = (e_(k) - e_(k-1)) / (s_(k) - s_(k-1))
    s_(k+1) = s_(k) - learning_rate*sigmoid(e_gradient)
    The energy at a certain scaling is calculated by using the simulation_function.

    Args:
        atoms(ASE atoms object): The configuration to find an optimal lattice constant for
        simulation_settings(dict): Dictionary of simulation settings to use.
        output_dict(dict): Dictionary to append the result to.
        learning_rate(float): The scaling converges quicker for larger values but if
             it's too large the scaling will oscillate and not converge at all

    Returns:
        scaling(float): Return the scaling factor which would give the inputed atoms object
            the lowest possible energy
    """
    old_scaling = 1
    scaling = 1.001
    e_scaling_gradient = 0
    number_of_iterations = 0
    old_energy_per_atom = atoms.get_total_energy()/len(atoms)
    best_scaling = 1
    best_energy_per_atom = old_energy_per_atom
    # Improve the lattice constant performing gradient descent until the gradient becomes
    # sufficiently small. A maximum of 100
    while (abs(e_scaling_gradient) > 0.01) or (number_of_iterations < 3):
        atoms_scaled = atoms.copy()
        # Perhaps we should calculate time average of energy instead of energy at the end
        # of each simulation.
        atoms_scaled.calc = atoms.calc
        atoms_scaled.set_cell(atoms.cell*scaling, scale_atoms=True)
        atoms_scaled, energy_per_atom = simple_simulation.run_simple_md_simulation(atoms_scaled, simulation_settings, True)
        e_scaling_gradient = (energy_per_atom-old_energy_per_atom)/(scaling-old_scaling)
        old_scaling = scaling
        scaling = old_scaling - learning_rate*e_scaling_gradient/(2+abs(e_scaling_gradient))
        print("Scaling gradient: ", e_scaling_gradient)
        print("Energy per atom: ", energy_per_atom)
        print("New scaling: ", scaling)
        old_energy_per_atom = energy_per_atom
        number_of_iterations += 1
        if energy_per_atom < best_energy_per_atom:
            best_scaling = old_scaling
            best_energy_per_atom = energy_per_atom
        if (number_of_iterations > 100):
            scaling = best_scaling
            break
    output_dict['optimal_scaling'].append(scaling)
    output_dict['iterations_to_find_scaling'].append(number_of_iterations)
    return scaling


if __name__ == "__main__":
    atoms = FaceCenteredCubic(size=(2, 2, 3), symbol='Cu', pbc=(1, 1, 0))
    atoms.calc = EMT()
    optimize_scaling(atoms, {'optimal_scaling': [], 'iterations_to_find_scaling': []})
