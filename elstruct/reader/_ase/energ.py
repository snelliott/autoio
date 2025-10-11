""" electronic energy readers
"""

from elstruct.par import Program
from elstruct.par import Method
from elstruct.par import Model
from elstruct.par import program_methods
from elstruct.par import program_models
from elstruct.par import method_is_mlip
from elstruct.par import mlip_from_method


PROG = Program.ASE

DOUB_HYB_DFT = [
    'b2plypd3'
]


def _hf_energy(output_dct):
    """ Reads the Hartree-Fock energy from the output dictionary

        :param output_dct: a dictionary of output information
        :type output_dct: dict
        :rtype: float
    """
    return output_dct['energy']


def _mp2_energy(output_dct):
    """ Reads the MP2 energy from the output file dictionary.
        Returns the energy in Hartrees.

        :param output_dct: dictionary of the program's output file
        :type output_dct: dict
        :rtype: float
    """
    return output_dct['energy']


def _dft_energy(output_dct):
    """ Reads the energy from most density functional theory methods
        from the output file dictionary. Returns the energy in Hartrees.

        :param output_dct: dictionary of the program's output file
        :type output_dct: dict
        :rtype: float
    """
    return output_dct['energy']


def _doub_hyb_dft_energy(output_dct):
    """ Read the energy from double-hybdrid density functional theory methods
        from the output file dictionary. Return the energy in Hartrees.

        :param output_dct: dictionary of the program's output file
        :type output_dct: dict
        :rtype: float
    """
    return output_dct['energy']


def _mlip_energy(output_dct):
    """ Reads the energy from machine learning interatomic potentials
        from the output file dictionary. Returns the energy in Hartrees.

        :param output_dct: dictionary of the program's output file
        :type output_dct: dict
        :rtype: float
    """
    return output_dct['energy']


# A dictionary of functions for reading the energy from the output, by method
ENERGY_READER_DCT = {
    (Method.HF[0], frozenset({})): _hf_energy,
    (Method.Corr.MP2[0], frozenset({})): _mp2_energy,
}
METHODS = program_methods(PROG)
MODELS = program_models(PROG)

# Add DFT methods to the reader dictionary
for METHOD in METHODS:
    if Method.is_standard_dft(METHOD):
        if METHOD not in DOUB_HYB_DFT:
            ENERGY_READER_DCT[(METHOD, frozenset({}))] = _dft_energy
        else:
            ENERGY_READER_DCT[(METHOD, frozenset({}))] = _doub_hyb_dft_energy
    elif Method.is_semi_empirical(METHOD):
        ENERGY_READER_DCT[(METHOD, frozenset({}))] = _dft_energy

for MODEL in MODELS:
    if Model.is_pretrained_model(MODEL):
        ENERGY_READER_DCT[(MODEL, frozenset({}))] = _mlip_energy
    elif Model.is_local_model(MODEL):
        ENERGY_READER_DCT[(MODEL, frozenset({}))] = _mlip_energy
# Check if we have added any unsupported methods to the energy reader
READ_METHODS = set(method[0] for method in ENERGY_READER_DCT)
assert READ_METHODS <= set(METHODS + MODELS)


def method_list():
    """ Constructs a list of available electronic structure methods.
    """
    return tuple(sorted(ENERGY_READER_DCT.keys()))


def energy(method, output_dct):
    """ Reads the the total electronic energy from the output file dictionary.
        Returns the energy in Hartrees.

        :param method: electronic structure method to read
        :type method: str
        :param output_dct: dictionary of the program's output file
        :type output_dct: dict 
        :rtype: float
    """
    # Parse the method and lists
    if method_is_mlip(method):
        core_method = mlip_from_method(method)
        pfxs = ()
    else:
        core_method, pfxs = Method.evaluate_method_type(method)
    full_method = (core_method, frozenset(pfxs))
    assert full_method in method_list()

    # Get the appropriate reader and call it
    if Model.contains(core_method):
        energy_reader = _mlip_energy
    elif Method.is_nonstandard_dft(core_method):
        if core_method not in DOUB_HYB_DFT:
            energy_reader = _dft_energy
        else:
            energy_reader = _doub_hyb_dft_energy
    else:
        energy_reader = ENERGY_READER_DCT[full_method]

    return energy_reader(output_dct)
