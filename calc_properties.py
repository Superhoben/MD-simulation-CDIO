from ase import units

def calc_temp(atoms):
    """ Calculates temperature of atoms object
    
    Args:
        atoms(ase atom object): the system to calculate the temperature for

    Returns:
        (float): the calculated temperature
    """
    ekin = atoms.get_kinetic_energy() / len(atoms)
    return ekin / (1.5 * units.kB)
