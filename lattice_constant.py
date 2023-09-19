from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT

def optimize_lattice_const(atoms):
    """ Find optimal scaling factor for the lattice constant given in the atoms object, 
        only searches through +-5%range with a 0.01% resolution

        Args: 
            atoms(ase atoms object): The configuration to find an optimal lattice constant for

        Returns: 
            lattice(float): Optimal lattice constant for the given configuration
    """

    print(atoms.cell)

    best_lattice_const = 0
    best_e_tot = 1000000
    atoms.calc = EMT()

    print(atoms.cell)
    for i in range(10):
        # +- 5% of initial size
        lattice_const = 0.01*i+0.95
        atoms1 = atoms.copy()
        atoms1.cell = atoms1.cell*lattice_const
        atoms1.calc = EMT()
        dyn = VelocityVerlet(atoms1, 5 * units.fs)
        dyn.run(200)
        epot = atoms1.get_potential_energy() / len(atoms1)
        ekin = atoms1.get_kinetic_energy() / len(atoms1)
        print("Etot = ")
        print(epot+ekin)
        

    #lattice_const

    #atoms1 = atoms
    #atoms1.cell = atoms1.cell*lattice_const

    print(atoms.cell)
    return




if __name__ == "__main__":
    optimize_lattice_const(FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [1,1,1]],
                          size=(2,2,3), symbol='Cu', pbc=(1,1,0)))