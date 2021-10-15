"""
Write various parts of a Chemkin mechanism file
"""

from chemkin_io.writer import reaction as writer_reac
from chemkin_io.writer import thermo as writer_therm
from chemkin_io.writer import _util as util


def write_chemkin_file(elem_tuple=None, mech_spc_dct=None, spc_nasa7_dct=None,
                       rxn_param_dct=None):
    """ Writes a Chemkin-formatted mechanism and/or thermo file. Writes
        the output to a text file.

        :param elem_tuple: tuple containing the element names
        :type elem_tuple: tuple
        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: {spc_name:data}
        :param spc_nasa7_dct: containing NASA-7 thermo data for each species
        :type spc_nasa7_dct: {spc_name:NASA-7 parameters}
        :param rxn_param_dct: containing the reaction parameters
        :type rxn_param_dct: {rxn:params}
    """

    total_str = ''
    if elem_tuple:
        elem_str = elements_block(elem_tuple)
        total_str += elem_str
    if mech_spc_dct:
        spc_str = species_block(mech_spc_dct)
        total_str += spc_str
    if spc_nasa7_dct:
        thermo_str = thermo_block(spc_nasa7_dct)
        total_str += thermo_str
    if rxn_param_dct:
        rxn_str = reactions_block(rxn_param_dct)
        total_str += rxn_str

    return total_str


def elements_block(elem_tuple):
    """ Writes the elements block of the mechanism file

        :param elem_tuple: tuple containing the element names
        :type elem_tuple: tuple
        :return elem_str: str containing the elements block
        :rtype: str
    """

    elem_str = 'ELEMENTS \n\n'
    for elem in elem_tuple:
        elem_str += elem + '\n'
    elem_str += '\nEND \n\n\n'

    return elem_str


def species_block(mech_spc_dct):
    """ Writes the species block of the mechanism file

        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: dct {spc_name:data}
        :return spc_str: str containing the species block
        :rtype: str
    """

    # Get the max species name length
    max_spc_len = 0
    for spc in mech_spc_dct.keys():
        if len(spc) > max_spc_len:
            max_spc_len = len(spc)

    # Get the max SMILES name length
    max_smiles_len = 0
    for spc, ident_dct in mech_spc_dct.items():
        if len(ident_dct['smiles']) > max_smiles_len:
            max_smiles_len = len(ident_dct['smiles'])

    buffer = 5

    # Write the spc_str
    spc_str = 'SPECIES \n\n'
    for spc, ident_dct in mech_spc_dct.items():
        spc_str += (
            '{0:<'+str(max_spc_len+buffer)+'s}{1:>9s}{2:<' +
            str(max_smiles_len+buffer)+'s}{3:>9s}{4:<9s}\n').format(
                spc, '! SMILES: ',
                ident_dct['smiles'], 'InChi: ', ident_dct['inchi'])

    spc_str += '\nEND \n\n\n'

    return spc_str


def thermo_block(spc_nasa7_dct):
    """ Writes the thermo block of the mechanism file

    """

    thermo_str = 'THERMO \n'
    thermo_str += '200.00    1000.00   5000.000  \n\n'
    for spc_name, params in spc_nasa7_dct.items():
        thermo_str += writer_therm.thermo_entry(spc_name, params)

    thermo_str += '\nEND\n\n\n'

    return thermo_str


def reactions_block(rxn_param_dct):
    """ Writes the reaction block of the mechanism file

        :param rxn_param_dct: dct containing the reaction parameters
        :type rxn_param_dct: dct {rxn: params}
        :return total_rxn_str: str containing the reaction block
        :rtype: str
    """

    # Get the length of the longest reaction name
    max_len = 0
    for rxn in rxn_param_dct.keys():
        rxn_name = util.format_rxn_name(rxn)
        if len(rxn_name) > max_len:
            max_len = len(rxn_name)

    rxn_str = 'REACTIONS     CAL/MOLE     MOLES\n\n'
    for rxn, params in rxn_param_dct.items():
        sing_rxn_str = writer_reac.get_ckin_str(rxn, params, max_len=max_len)
        rxn_str += sing_rxn_str
    rxn_str += '\n\nEND\n\n'

    return rxn_str
