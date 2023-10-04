from ase.lattice.cubic import FaceCenteredCubic
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase import Atoms


def optimize_lattice_const(atoms):
    """ Finds optimal scaling factor for the lattice constant given in the atoms object,
        only searches through +-5% range with a 0.002% resolution.

        Args:
            atoms(ase atoms object): The configuration to find an optimal lattice constant for

        Returns:
            best_lattice_const_scaling(float): The optimal scaling for the cell of the atoms object
    """
    # Initialize values
    best_lattice_const_scaling = 0
    best_e_tot = float('inf')
    original_cell = atoms.cell
    # Maybe add ability to choose the potential?
    atoms.calc = EMT()

    for i in range(101):
        # +- 5% of initial size, this is chosen arbitrarily so maybe something else is better?
        lattice_const = 0.01 * i + 0.95
        atoms_copy = atoms.copy()
        atoms_copy.calc = EMT()
        atoms_copy.set_cell(original_cell*lattice_const, scale_atoms=True)
        # No evolution in time is needed since energy is preserved in a closed system, however this
        # is keept as reminder for when we use a diffrent ensemble where energy can be transferred.
        dyn = VelocityVerlet(atoms_copy, 5 * units.fs)
        dyn.run(1)
        e_tot = atoms_copy.get_total_energy()
        if e_tot < best_e_tot:
            best_e_tot = e_tot
            best_lattice_const_scaling = lattice_const
    return best_lattice_const_scaling


def example_simulation_function(atoms: Atoms):
    """An example of a function which evolves an atoms object which can
    be used by the optimize_lattice_const_gradient_descent function.

    Args:
        atoms(ASE atoms object): The intial configuraion to evolve

    Returns:
        atoms(ASE atoms object): The final evolved version
    """
    atoms.calc = EMT()
    dyn = VelocityVerlet(atoms, 5 * units.fs)
    # No evolution in time is needed since energy is preserved in a closed system, however this
    # is keept as reminder for when we use a diffrent ensemble where energy can be transferred.
    dyn.run(1)
    return atoms


def optimize_lattice_const_gradient_descent(atoms, simulation_function, learning_rate=0.05):
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
        learning_rate(float): The scaling converges quicker for larger values but if
             it's too large the scaling will oscillate and not converge at all
        simulation_function(function[atoms_object]->atoms_object): A function which take an
            atoms obejct, evolves this object throughout time and return the atoms object

    Returns:
        scaling(float): Return the scaling factor which would give the inputed atoms object
            the lowest possible energy
        energy(float): Also return the energy for the scaled atoms object
        number_of_iterations(int): Shows how many iterations it took for the energy to converge
    """
    old_scaling = 1
    scaling = 1.001
    e_scaling_gradient = 0
    number_of_iterations = 0
    old_energy_per_atom = simulation_function(atoms.copy()).get_total_energy()
    # Improve the lattice constant performing gradient descent until the gradient becomes
    # sufficiently small. At least 3 iteration will be performed even if the gradient
    # start small.
    while (abs(e_scaling_gradient) > 0.01) or (number_of_iterations < 3):
        atoms_scaled = atoms.copy()
        atoms_scaled.set_cell(atoms.cell*scaling, scale_atoms=True)
        energy_per_atom = simulation_function(atoms_scaled).get_total_energy()/len(atoms_scaled)
        e_scaling_gradient = (energy_per_atom-old_energy_per_atom)/(scaling-old_scaling)
        old_scaling = scaling
        scaling = old_scaling - learning_rate*e_scaling_gradient/(2+abs(e_scaling_gradient))
#        print("Scaling gradient: ", e_scaling_gradient)
#        print("Energy per atom: ", energy_per_atom)
#        print("New scaling: ", scaling)
        old_energy_per_atom = energy_per_atom
        number_of_iterations += 1
    return scaling, energy_per_atom, number_of_iterations


if __name__ == "__main__":
    optimize_lattice_const(FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]],
                                             size=(2, 2, 3), symbol='Cu', pbc=(1, 1, 0)))
