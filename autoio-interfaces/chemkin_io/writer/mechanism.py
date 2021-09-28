"""
Write various parts of a Chemkin mechanism file
"""

from chemkin_io.writer import reaction as writer_reac
from chemkin_io.writer import thermo as writer_therm
from chemkin_io.writer import _util as util


def write_chemkin_file(elem_tuple=None, spc_dct=None, spc_nasa7_dct=None,
                       rxn_param_dct=None, comments=None):
    """ Writes a Chemkin-formatted mechanism and/or thermo file. Writes
        the output to a text file.

        :param elem_tuple: tuple containing the element names
        :type elem_tuple: tuple
        :param spc_dct: data
        :type spc_dct: {spc_name:data}
        :param spc_nasa7_dct: containing NASA-7 thermo data for each species
        :type spc_nasa7_dct: {spc_name:NASA-7 parameters}
        :param rxn_param_dct: containing the reaction parameters
        :type rxn_param_dct: {rxn:params}
    """

    total_str = ''
    if elem_tuple:
        elem_str = elements_block(elem_tuple)
        total_str += elem_str
    if spc_dct:
        spc_str = species_block(spc_dct)
        total_str += spc_str
    if spc_nasa7_dct:
        thermo_str = thermo_block(spc_nasa7_dct)
        total_str += thermo_str
    if rxn_param_dct:
        rxn_str = reactions_block(rxn_param_dct, comments=comments)
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


def species_block(spc_ident_dct):
    """ Writes the species block of the mechanism file

        :param spc_dct: dct containing the species data
        :type spc_dct: dct {spc_name:data}
        :return spc_str: str containing the species block
        :rtype: str
    """
    # Get the max species name length
    max_spc_len = 0
    for spc in spc_ident_dct.keys():
        if len(spc) > max_spc_len:
            max_spc_len = len(spc)

    # Get the max SMILES name length
    max_smiles_len = 0
    for spc, ident_dct in spc_ident_dct.items():
        if len(ident_dct['smiles']) > max_smiles_len:
            max_smiles_len = len(ident_dct['smiles'])

    buffer = 5

    # Write the spc_str
    spc_str = 'SPECIES \n\n'
    for spc, ident_dct in spc_ident_dct.items():
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


def reactions_block(rxn_param_dct, comments=None):

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


def reactions_block_old(rxn_param_dct, comments=None):
    """ Writes the reaction block of the mechanism file

        :param rxn_param_dct: dct containing the reaction parameters
        :type rxn_param_dct: dct {rxn:params}
        :return total_rxn_str: str containing the reaction block
        :rtype: str
    """

    # Get the length of the longest reaction name
    max_len = 0
    for rxn, param_dct in rxn_param_dct.items():
        rxn_name = util.format_rxn_name(rxn)
        if len(rxn_name) > max_len:
            max_len = len(rxn_name)

    # Loop through each reaction and get the string to write to text file
    total_rxn_str = 'REACTIONS     CAL/MOLE     MOLES\n\n'
    for rxn, param_dct in rxn_param_dct.items():

        # preprocess potential duplicate PLOGs written as separate tuples:
        param_dct = util.merge_plog_dct(param_dct)
        # loop over the parameter entries:
        # might be a single one or a duplicate one
        for param_set_i, param_dct_vals in enumerate(param_dct):
            # Convert the reaction name from tuple of tuples to string
            # (Note: this includes '+M' or '(+M)' if appropriate)
            rxn_name = util.format_rxn_name(rxn)

            if param_dct_vals[3] is not None:  # Chebyshev
                assert param_dct_vals[0] is not None, (
                    f'For {rxn}, Chebyshev params included but' +
                    ' highP params absent'
                )
                # this spot is usually high-P params,
                # but is instead 1-atm for Chebyshev
                one_atm_params = param_dct_vals[0]
                alpha = param_dct_vals[3]['alpha_elm']
                tmin = param_dct_vals[3]['t_limits'][0]
                tmax = param_dct_vals[3]['t_limits'][1]
                pmin = param_dct_vals[3]['p_limits'][0]
                pmax = param_dct_vals[3]['p_limits'][1]
                rxn_str = writer_reac.chebyshev(
                    rxn_name, one_atm_params, alpha,
                    tmin, tmax, pmin, pmax, max_length=max_len)

            elif param_dct_vals[4] is not None:  # PLOG
                plog_dct = param_dct_vals[4]
                rxn_str = writer_reac.plog(
                    rxn_name, plog_dct, max_length=max_len)

            elif param_dct_vals[2] is not None:  # Troe
                assert param_dct_vals[0] is not None, (
                    f'For {rxn}, Troe params included, highP params absent'
                )
                assert param_dct_vals[1] is not None, (
                    f'For {rxn}, Troe, highP params included,' +
                    ' lowP params absent'
                )
                assert rxn[2][0] is not None, (
                    f'For {rxn}, Troe, highP, lowP params included,' +
                    ' third body absent'
                )

                highp_params = param_dct_vals[0]
                lowp_params = param_dct_vals[1]
                troe_params = param_dct_vals[2]
                collid_factors = param_dct_vals[5]
                rxn_str = writer_reac.troe(
                    rxn_name, highp_params, lowp_params, troe_params,
                    colliders=collid_factors, max_length=max_len
                )

            elif param_dct_vals[1] is not None:  # Lindemann
                assert param_dct_vals[0] is not None, (
                    f'For {rxn}, lowP params included, highP params absent'
                )
                assert rxn[2][0] is not None, (
                    f'For {rxn}, highP, lowP params included,' +
                    ' third body absent'
                )
                highp_params = param_dct_vals[0]
                lowp_params = param_dct_vals[1]
                collid_factors = param_dct_vals[5]
                rxn_str = writer_reac.lindemann(
                    rxn_name, highp_params, lowp_params,
                    colliders=collid_factors,
                    max_length=max_len
                )

            else:  # Simple Arrhenius
                assert param_dct_vals[0] is not None, (
                    f'For {rxn}, the highP params absent'
                )
                highp_params = param_dct_vals[0]
                collid_factors = param_dct_vals[5]
                rxn_str = writer_reac.arrhenius(
                    rxn_name, highp_params,
                    colliders=collid_factors, max_length=max_len)

            if comments and param_set_i == 0:
                # add inline comments on the first line
                if isinstance(comments[rxn], dict):
                    if comments[rxn]['cmts_inline'] != '':
                        rxn_str_split = rxn_str.split('\n')
                        rxn_str_split[0] = rxn_str_split[0] + \
                            ' ' + comments[rxn]['cmts_inline']
                        # rewrite rxn_str
                        rxn_str = '\n'.join(rxn_str_split)

                # check for comments: header
                if isinstance(comments[rxn], dict):
                    total_rxn_str += comments[rxn]['cmts_top']

            if len(param_dct) > 1:
                # duplicate reactions
                rxn_str += 'DUP \n'

            total_rxn_str += rxn_str

    total_rxn_str += '\nEND \n'

    return total_rxn_str
