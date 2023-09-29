from ase.lattice.cubic import FaceCenteredCubic
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT


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


if __name__ == "__main__":
    optimize_lattice_const(FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]],
                                             size=(2, 2, 3), symbol='Cu', pbc=(1, 1, 0)))
