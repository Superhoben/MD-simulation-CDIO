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
    atoms1 = atoms.copy()
    # Maybe add ability to choose the potential?
    atoms1.calc = EMT()
    original_cell = atoms.cell
    dyn = VelocityVerlet(atoms1, 5 * units.fs)

    for i in range(51):
        # +- 5% of initial size, this is chosen arbitrarily so maybe something else is better?
        lattice_const = 0.002 * i + 0.95
        # Perhaps it's better to not "reset" the cell and instead do atoms1.cell = atoms1.cell + 0.002*original_cell?
        # Then of course row above must be removed and another value for atoms1.cell must be initialized before loop.
        atoms1.cell = original_cell * lattice_const
        # We chose this arbitrarily, perhaps other number is better?
        dyn.run(400)
        # Perhaps unnecessary to divide with length? However, if big cell then might lose accuracy in comparison
        # Perhaps this should be sampled over time also?
        e_tot = (atoms1.get_potential_energy() + atoms1.get_kinetic_energy()) / len(atoms1)
        if e_tot < best_e_tot:
            best_e_tot = e_tot
            best_lattice_const_scaling = lattice_const

    return best_lattice_const_scaling


def example_simulation_function(atoms):
    """An example of a function which evolves an atoms object which can
    be used by the optimize_lattice_const_gradient_descent function.

    Args:
        atoms(ASE atoms object): The intial configuraion to evolve

    Returns:
        atoms(ASE atoms object): The final evolved version
    """
    atoms.calc = EMT()
    dyn = VelocityVerlet(atoms, 5 * units.fs)
    dyn.run(100)
    return atoms


def optimize_lattice_const_gradient_descent(atoms, learning_rate, simulation_function):
    """Finds the optimal sclang of the lattice constant using the gradient descent method
    modified by a sigmoid function since energy changes can become very extrem as some
    potential depend on the -12th power of the distance so starting with half the correct
    distance can result in approximately 4000 times the energy which would make a normal
    gradient search go crazy even with a reasonable learning rate. The new scaling s_(k+1) 
    is calculated with using the change in energy, e, as:
    e_gradient = (e_(k) - e_(k-1)) / (s_(k) - s_(k-1))
    s_(k+1) = s_(k) - learning_rate*sigmoid(e_gradient) 
    The energy at a certain scaling is calculated by using the simulation_function.

    Args:
        atoms(ASE atoms object): The configuration to find an optimal lattice constant for
        learning_rate (_type_): The scaling converges quicker for larger values but if
             it's too large the scaling will oscillate and not converge at all
        simulation_function(function[atoms_object]->atoms_object): A function which take an
            atoms obejct, evolves this object throughout time and return the atoms object


    Returns:
        scaling(float): Return the scaling factor which would give the inputed atoms object
            the lowest possible energy 
        energy(float): Also return the energy for the scaled atoms object
        number_of_iterations(int): Shows how many iterations it took for the energy to converge
    """
    old_energy_per_atom = simulation_function(atoms)
    old_scaling = 1
    scaling = 1.1
    e_scaling_gradient = 0
    number_of_iterations = 0
    while (e_scaling_gradient < 0.01) or (number_of_iterations < 3):
        atoms_scaled = atoms.copy()
        atoms_scaled.cell = atoms.cell*scaling
        energy_per_atom = simulation_function(atoms_scaled).get_total_energy()/len(atoms)
        e_scaling_gradient = (energy_per_atom-old_energy_per_atom)/(scaling-old_scaling)
        scaling = scaling - learning_rate*e_scaling_gradient/(1+abs(e_scaling_gradient))
        old_energy_per_atom = energy_per_atom
        number_of_iterations += 1
    return scaling, energy_per_atom, number_of_iterations


if __name__ == "__main__":
    optimize_lattice_const(FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]],
                                             size=(2, 2, 3), symbol='Cu', pbc=(1, 1, 0)))
