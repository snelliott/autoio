"""Functions to run ASE calculators."""

from ase import Atoms
from ase import units
from ase.optimize import BFGS
# from sella import Sella

def get_calculator(calculator, family=None):
    """Get a new instance of the appropriate ASE calculator.
    
    Args:
        program: The program name ('psi4', 'mace_mp', etc.)
    Returns:
        The calculator class
    Raises:
        ValueError: If the program is not supported
    """
    # Lazy import the appropriate calculator
    if calculator == 'psi4':
        import psi4
        from ase.calculators.psi4 import Psi4
        return Psi4
    elif calculator == 'mace':
        if family == 'mace_mp':
            from mace.calculators import mace_mp
            return mace_mp
        elif family == 'mace_off':
            from mace.calculators import mace_off
            return mace_off
        elif family == 'mace_anicc':
            from mace.calculators import mace_anicc
            return mace_anicc
        else:
            raise ValueError(f"Unsupported MACE family: {family}")
    # elif calculator == 'uma':
    # elif calculator == 'nwx':
    #     from ase.calculators.nwchem import NWChem
    #     return NWChem


def from_calc_dictionary(input_dct, script_str):
    """ Run an ASE calculation from an input dictionary and return results

    Args:
        input_dct: Dictionary containing:
            - calculator: Name of ASE calculator to use ('psi4' or 'nwx')
            - atom_parameters:
                - symbols: List of atomic symbols
                - positions: List of [x,y,z] coordinates
            - basis: Basis set name
            - method: Quantum chemistry method
            - charge: Molecular charge
            - multiplicity: Spin multiplicity
            - reference: Reference type ('rhf', 'uhf', etc.)
        script_str: String of script that contains ase_<calculator>
    Returns:
        dict: Results containing:
            - energy: Total energy in eV
            - forces: Forces on atoms in eV/Å
            - dipole: Dipole moment in e⋅Å
            - charges: Atomic charges
            - parameters: Input parameters used
            - version: Version string of the program used
    """
    # Set up Atoms object
    atom_parameters = input_dct.get('atom_parameters', {})
    job = input_dct.get('job', 'energy')
    atoms = Atoms(
        symbols=atom_parameters.get('symbols', []),
        positions=atom_parameters.get('positions', []),
    )
    # Set up calculator
    calculator = None
    if len(script_str.split('ase_')) > 1:
        calculator = script_str.split('ase_')[-1].strip()
        if 'method' in input_dct.keys():
            kwargs = {
                'atoms': atoms,
                'basis': input_dct.get('basis'),
                'method': input_dct.get('method'),
                'charge': input_dct.get('charge', 0),
                'multiplicity': input_dct.get('multiplicity', 1),
                'reference': input_dct.get('reference')
            }
            family = None
        elif 'family' in input_dct.keys():
            family = input_dct.get('family')
            if 'mlip' in input_dct.keys():
                kwargs = {
                    'model': input_dct.get('mlip'),
                    'device': 'cpu' 
                }
            else:
                kwargs = {}
    calc = get_calculator(calculator, family)(**kwargs)
    
    # Attach calculator to atoms and run
    atoms.calc = calc
    
    # Run calculation and get basic properties that all calculators support
    if job == 'energy':
        energy = atoms.get_potential_energy() / units.Hartree
        positions = atoms.get_positions()
    elif job == 'optimize':
        if not input_dct.get('saddle', False):
            dyn = BFGS(atoms)
        # else:
        #     dyn = Sella(atoms)
        dyn.run(
            fmax=input_dct.get('gconv', 0.05),
            steps=input_dct.get('maxcyc', 100),
            opt_type=input_dct.get('opt_type', 'standard')
            ) 
        energy = atoms.get_potential_energy()
        positions = atoms.get_positions()

    # Initialize results with guaranteed properties
    results = {
        'energy': energy,
        'positions': positions,
        # 'forces': forces.tolist(),
        'parameters': dict(calc.parameters)  # Convert to regular dict for JSON serialization
    }

    # Add calculator-specific properties based on calculator type
    if calculator == 'psi4':
        results['version'] = f'ase-psi4'
        results['normal_termination'] = True
        # Psi4 supports both dipole moments and charges
        # if 'properties' not in calc.parameters or 'dipole' in calc.parameters.get('properties', []):
        #     results['dipole'] = atoms.get_dipole_moment().tolist()
        # if 'properties' not in calc.parameters or 'mulliken' in calc.parameters.get('properties', []):
        #     results['charges'] = atoms.get_charges().tolist()
    elif 'mace' in calculator:
        # MACE specific properties
        results['version'] = f'ase-mace'
        results['normal_termination'] = True
        # MACE does not provide dipole moments or charges
        # if 'properties' not in calc.parameters or 'dipole' in calc.parameters.get('properties', []):
        #     results['dipole'] = atoms.get_dipole_moment().tolist()
        # if 'properties' not in calc.parameters or 'mulliken' in calc.parameters.get('properties', []):
        #     results['charges'] = atoms.get_charges().tolist()
    elif calculator == 'nwx':
        # NWChem specific properties
        pass
        # if 'dipole_moment' in calc.get_implemented_properties():
        #     results['dipole'] = atoms.get_dipole_moment().tolist()
        # if 'charges' in calc.get_implemented_properties():
        #     results['charges'] = atoms.get_charges().tolist()

    return results
    