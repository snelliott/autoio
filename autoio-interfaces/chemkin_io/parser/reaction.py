""" functions operating on the reactions block string
"""

import sys
import collections
import itertools
import math
import copy
import numpy as np
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from ioformat import headlined_sections
from phydat import phycon


# Various strings needed to parse the data sections of the Reaction block
CHEMKIN_ARROW = (app.maybe(app.escape('<')) + app.escape('=') +
                 app.maybe(app.escape('>')))
CHEMKIN_PLUS_EM = app.PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = app.escape('(') + app.PLUS + 'M' + app.escape(')')
CHEMKIN_PAREN_PLUS = app.escape('(') + app.PLUS
CHEMKIN_PAREN_CLOSE = app.escape(')')

SPECIES_NAME_PATTERN = (
    r'[^\s=+\-]' +
    app.zero_or_more(app.one_of_these(
        [app.LETTER, app.DIGIT, r'[#,()\-_]',
         app.escape('*'), app.escape('(+)'),
         app.escape('['), app.escape(']')])) +
    app.zero_or_more(app.PLUS)
)
SPECIES_NAMES_PATTERN = app.series(
    app.padded(SPECIES_NAME_PATTERN), app.padded(app.PLUS))

REACTION_PATTERN = (SPECIES_NAMES_PATTERN + app.padded(CHEMKIN_ARROW) +
                    SPECIES_NAMES_PATTERN)
COEFF_PATTERN = (app.NUMBER + app.LINESPACES + app.NUMBER +
                 app.LINESPACES + app.NUMBER)
COMMENTS_PATTERN = app.escape(
    '!') + app.capturing(app.one_or_more(app.WILDCARD2))

BAD_STRS = ['inf', 'INF', 'nan']


def param_dct(block_str, ea_units, a_units):
    """ Parses all of the chemical equations and corresponding fitting
        parameters in the reactions block of the mechanism input file
        and subsequently pulls all of the species names and fitting
        parameters from the data string; this information is stored in a list.

        :param block_str: string for reactions block
        :type block_str: str
        :return rxn_param_dct: dict{(reacs,prods): param_list}
        :rtype: dict
    """

    rxn_dstr_lst = data_strings(block_str)

    if rxn_dstr_lst is not None:

        # Create an iterator that repeats the units inputs
        many_ea_units = list(
            itertools.repeat(ea_units, times=len(rxn_dstr_lst)))
        many_a_units = list(
            itertools.repeat(a_units, times=len(rxn_dstr_lst)))

        reac_and_prods = list(zip(
            map(reactant_names, rxn_dstr_lst),
            map(product_names, rxn_dstr_lst),
            map(third_body, rxn_dstr_lst)))

        params = list(zip(
            map(high_p_parameters, rxn_dstr_lst, many_ea_units, many_a_units),
            map(low_p_parameters, rxn_dstr_lst, many_ea_units, many_a_units),
            map(troe_parameters, rxn_dstr_lst),
            map(chebyshev_parameters, rxn_dstr_lst, many_a_units),
            map(plog_parameters, rxn_dstr_lst, many_ea_units, many_a_units),
            map(collider_enhance_factors, rxn_dstr_lst)))

        # Fix any duplicates in Arrhenius or PLOG reactions
        reac_and_prods, params = fix_duplicates(reac_and_prods, params)
        rxn_param_dct = dict(zip(reac_and_prods, params))

    else:
        rxn_param_dct = None

    return rxn_param_dct


def data_strings(block_str, remove_bad_fits=False):
    """ Parses all of the chemical equations and corresponding fitting
        parameters in the reactions block of the mechanism input file
        and stores them in a list.

        :param block_str: string for reactions block
        :type block_str: str
        :param remove_bad_fits: remove reactions with bad fits
        :type remove_bad_fits: bool
        :return rxn_dstrs: strings containing eqns and params for all reactions
        :rtype: list(str)
    """

    if block_str.strip():
        rxn_dstrs = headlined_sections(
            string=block_str.strip(),
            headline_pattern=CHEMKIN_ARROW
        )
        if remove_bad_fits:
            rxn_dstrs = [dstr for dstr in rxn_dstrs
                         if not any(string in dstr for string in BAD_STRS)]
    else:
        rxn_dstrs = None

    return rxn_dstrs


def reactant_names(rxn_dstr):
    """ Parses the data string for a reaction in the reactions block
        for the line containing the chemical equation in order to
        read the names of the reactant species.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return names: names of the reactants
        :rtype: tuple(str)
    """
    pattern = _first_line_pattern(
        rct_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        prd_ptt=SPECIES_NAMES_PATTERN,
        param_ptt=COEFF_PATTERN
    )
    string = apf.first_capture(pattern, rxn_dstr)
    try:
        names = _split_reagent_string(string)
    except TypeError:
        print('Reaction line not formatted correctly:\n', rxn_dstr)
        print('Check that there are three numbers after the rxn equation')
        sys.exit()

    return names


def product_names(rxn_dstr):
    """ Parses the data string for a reaction in the reactions block
        for the line containing the chemical equation in order to
        read the names of the product species.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return names: names of the products
        :rtype: tuple(str)
    """
    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        param_ptt=COEFF_PATTERN
    )
    string = apf.first_capture(pattern, rxn_dstr)
    try:
        names = _split_reagent_string(string)
    except TypeError:
        print('Reaction line not formatted correctly:\n', rxn_dstr)
        print('Check that there are three numbers after the rxn equation')
        sys.exit()

    return names


def third_body(rxn_dstr):
    """ Parses the data string for a reaction in the reactions block
        for the line containing the chemical equation in order to
        read the names of the third body collider if present

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return trd_body: names of the colliders and corresponding fraction
        :rtype: tuple(str)
    """
    pattern = _first_line_pattern(
        rct_ptt=app.capturing(SPECIES_NAMES_PATTERN),
        prd_ptt=SPECIES_NAMES_PATTERN,
        param_ptt=app.maybe(COEFF_PATTERN)
    )
    rgt_str = apf.first_capture(pattern, rxn_dstr)
    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_split_paren = apf.split(CHEMKIN_PAREN_PLUS, rgt_str)
    rgt_split_plus = apf.split(app.PLUS, rgt_str)

    if len(rgt_split_paren) > 1:
        trd_body = '(+' + apf.split(CHEMKIN_PAREN_CLOSE,
                                    rgt_split_paren[1])[0] + ')'

    elif 'M' in rgt_split_plus:
        trd_body = '+M'

    else:
        trd_body = None

    trd_body = (trd_body,)

    return trd_body


def high_p_parameters(rxn_dstr, ea_units, a_units):
    """ Parses the data string for a reaction in the reactions block
        for the line containing the chemical equation in order to
        read the fitting parameters that are on the same line.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: Arrhenius fitting parameters for high-P rates
        :rtype: list(float)
    """
    pattern = _first_line_pattern(
        rct_ptt=SPECIES_NAMES_PATTERN,
        prd_ptt=SPECIES_NAMES_PATTERN,
        param_ptt=app.capturing(COEFF_PATTERN)
    )

    string_lst = apf.all_captures(pattern, rxn_dstr)
    if string_lst:
        fake_params = []
        for string in string_lst:
            fake_params.append(list(ap_cast(string.split())))
            params = fake_params[0]

        # Convert the units of Ea and A
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_dstr, a_units)
        params[2] = params[2] * ea_conv_factor
        params[0] = params[0] * a_conv_factor

    else:
        params = None

    return params


def low_p_parameters(rxn_dstr, ea_units, a_units):
    """ Parses the data string for a reaction in the reactions block
        for a line containing the low-pressure fitting parameters,
        then reads the parameters from this line.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: Arrhenius fitting parameters for low-P rates
        :rtype: list(float)
    """
    pattern = (
        'LOW' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.zero_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    cap1 = apf.first_capture(pattern, rxn_dstr)
    if cap1 is not None:
        params = [float(val) for val in cap1]

        # Convert the units of Ea and A
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_dstr, a_units)
        params[2] = params[2] * ea_conv_factor
        params[0] = params[0] * a_conv_factor

    else:
        params = None

    return params


def troe_parameters(rxn_dstr):
    """ Parses the data string for a reaction in the reactions block
        for a line containing the Troe fitting parameters,
        then reads the parameters from this line.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: Troe fitting parameters
        :rtype: list(float)
    """
    pattern = (
        'TROE' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.zero_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.maybe(app.one_or_more(app.SPACE) + app.capturing(app.NUMBER)) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    cap1 = apf.first_capture(pattern, rxn_dstr)

    if cap1 is not None:
        params = []
        for val in cap1:
            if val is not None:
                params.append(float(val))
            else:
                params.append(None)
    else:
        params = None

    return params


def chebyshev_parameters(rxn_dstr, a_units='moles'):
    """ Parses the data string for a reaction in the reactions block
        for the lines containing the Chebyshevs fitting parameters,
        then reads the parameters from these lines.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: Chebyshev fitting parameters
        :rtype: dict[param: value]
    """
    original_rxn_dstr = rxn_dstr
    rxn_dstr = apf.remove(COMMENTS_PATTERN, rxn_dstr)

    tcheb_pattern = (
        'TCHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.zero_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    pcheb_pattern = (
        'PCHEB' + app.zero_or_more(app.SPACE) + app.escape('/') +
        app.zero_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    cheb_pattern = (
        app.not_preceded_by(app.one_of_these(['T', 'P'])) +
        'CHEB' + app.zero_or_more(app.SPACE) +
        app.escape('/') + app.capturing(app.one_or_more(app.WILDCARD2)
                                        ) + app.escape('/')
    )

    cheb_params_raw = apf.all_captures(cheb_pattern, rxn_dstr)

    if cheb_params_raw:
        params = {}
        # Get the temp and pressure limits;
        # add the Chemkin default values if they don't exist
        cheb_temps = apf.first_capture(tcheb_pattern, rxn_dstr)
        cheb_pressures = apf.first_capture(pcheb_pattern, rxn_dstr)
        if cheb_temps is None:
            cheb_temps = ('300.00', '2500.00')
            print(
                'No Chebyshev temperature limits specified' +
                ' for the below reaction.' +
                f' Assuming 300 and 2500 K. \n \n {original_rxn_dstr}\n')
        if cheb_pressures is None:
            cheb_pressures = ('0.001', '100.00')
            print(
                'No Chebyshev pressure limits specified' +
                ' for the below reaction.' +
                f' Assuming 0.001 and 100 atm. \n \n {original_rxn_dstr}\n')

        # Get all the numbers from the CHEB parameters
        cheb_params = []
        for cheb_line in cheb_params_raw:
            cheb_params.extend(cheb_line.split())

        # Get the cheb array dimensions N and M, which are the first two
        # entries of the CHEB params
        cheb_n = int(math.floor(float(cheb_params[0])))
        cheb_m = int(math.floor(float(cheb_params[1])))

        # Start on the third value (after N and M)
        #  and get all the polynomial coefficients
        coeffs = []
        for idx, coeff in enumerate(cheb_params[2:]):
            # extra coefficients are allowed but ignored
            if idx+1 > (cheb_n*cheb_m):
                break
            coeffs.append(coeff)
        assert len(coeffs) == (cheb_n*cheb_m), (
            f'For the below reaction, there should be {cheb_n*cheb_m}' +
            ' Chebyshev polynomial' +
            f' coefficients, but there are only {len(coeffs)}.' +
            f' \n \n {original_rxn_dstr}\n'
        )
        alpha = np.array(list(map(float, coeffs)))

        params['t_limits'] = [float(val) for val in cheb_temps]
        params['p_limits'] = [float(val) for val in cheb_pressures]
        params['alpha_elm'] = alpha.reshape([cheb_n, cheb_m])
        params['a_units'] = a_units

    else:
        params = None

    return params


def plog_parameters(rxn_dstr, ea_units, a_units):
    """ Parses the data string for a reaction in the reactions block
        for the lines containing the PLOG fitting parameters,
        then reads the parameters from these lines.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: PLOG fitting parameters
        :rtype: dict[pressure: params]
    """
    pattern = (
        'PLOG' +
        app.zero_or_more(app.SPACE) + app.escape('/') +
        app.zero_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.one_or_more(app.SPACE) + app.capturing(app.NUMBER) +
        app.zero_or_more(app.SPACE) + app.escape('/')
    )
    params_lst = apf.all_captures(pattern, rxn_dstr)

    # Build dictionary of parameters, indexed by parameter
    if params_lst:

        # Get the Ea and A conversion factors
        ea_conv_factor = get_ea_conv_factor(ea_units)
        a_conv_factor = get_a_conv_factor(rxn_dstr, a_units)
        params = {}
        for param in params_lst:
            pressure = float(param[0])
            vals = list(map(float, param[1:]))
            vals[2] = vals[2] * ea_conv_factor
            vals[0] = vals[0] * a_conv_factor
            if pressure not in params:
                params[pressure] = vals
            else:
                params[pressure].extend(vals)  # add duplicate expressions

    else:
        params = None

    return params


def collider_enhance_factors(rxn_dstr):
    """ Parses the data string for a reaction in the reactions block
        for the line containing the names of several bath gases and
        their corresponding collision enhancement factors.

        :param rxn_dstr: data string for species in reaction block
        :type rxn_dstr: str
        :return params: Collision enhanncement factors for each bath gas
        :rtype: dict[bath name: enhancement factors]
    """
    bad_strings = ('DUP', 'LOW', 'TROE', 'CHEB', 'PLOG', CHEMKIN_ARROW)
    species_char = app.one_of_these([
        app.LETTER, app.DIGIT,
        app.escape('('), app.escape(')'),
        app.UNDERSCORE])
    species_name = app.one_or_more(species_char)

    # Loop over the lines and search for string with collider facts
    if ('LOW' in rxn_dstr or 'TROE' in rxn_dstr
            or 'M=' in rxn_dstr or 'M =' in rxn_dstr):
        params = {}
        for line in rxn_dstr.splitlines():
            if not any(apf.has_match(string, line) for string in bad_strings):
                factor_pattern = (
                    app.capturing(species_name) + app.zero_or_more(app.SPACE) +
                    app.escape('/') + app.zero_or_more(app.SPACE) +
                    app.capturing(app.NUMBER) + app.zero_or_more(app.SPACE) +
                    app.escape('/') + app.zero_or_more(app.SPACE)
                )
                baths = apf.all_captures(factor_pattern, line)
                if baths:
                    for bath in baths:
                        params[bath[0]] = float(bath[1])
        # If nothing was put into the dictionary, set it to None
        if not params:
            params = None
    else:
        params = None

    return params


def _first_line_pattern(rct_ptt, prd_ptt, param_ptt):
    """ Defines the pattern for the first line in a reaction data
        string that contains the chemical equation and high-pressure
        fitting parameters for the reaction.

        :param rct_ptt: string pattern for the reactant species
        :type rct_ptt: str
        :param prd_ptt: string pattern for the product species
        :type prd_ptt: str
        :param param_ptt: string pattern for high-pressure parameters
        :type param_ptt: str
        :rtype: str
    """
    return (rct_ptt + app.padded(CHEMKIN_ARROW) + prd_ptt +
            app.LINESPACES + param_ptt)


def _split_reagent_string(rgt_str):
    """ Parses out the names of all the species given in a string with
        the chemical equation within the reactions block.

        :param rgt_str: string with the reaction chemical equation
        :type rgt_str: str
        :return rgts: names of the species in the reaction
        :type rgts: tuple(str)
    """

    def _interpret_reagent_count(rgt_cnt_str):
        """ Count the species in a string containing one side
            of a chemical equation.

            :param rgt_cnt_str: string of one side of chemcial equation
            :type rgt_cnt_str: str
            :return: rgts: names of species from string
            :rtype: tuple(str)
        """
        _pattern = (app.STRING_START + app.capturing(app.maybe(app.DIGIT)) +
                    app.capturing(app.one_or_more(app.NONSPACE)))
        cnt, rgt = apf.first_capture(_pattern, rgt_cnt_str)
        cnt = int(cnt) if cnt else 1
        rgts = (rgt,) * cnt
        return rgts

    rgt_str = apf.remove(app.LINESPACES, rgt_str)
    rgt_str = apf.remove(CHEMKIN_PAREN_PLUS_EM, rgt_str)
    rgt_split = apf.split(CHEMKIN_PAREN_PLUS, rgt_str)

    if len(rgt_split) > 1:
        rgt_str = rgt_split[0]

    pattern = app.PLUS + app.not_followed_by(app.PLUS)
    rgt_cnt_strs = apf.split(pattern, rgt_str)

    rgts = tuple(itertools.chain(*map(_interpret_reagent_count, rgt_cnt_strs)))
    # remove '+M':
    if 'M' in rgts:
        rgts_lst = list(rgts)
        rgts_lst.remove('M')
        rgts = tuple(rgts_lst)

    return rgts


def fix_duplicates(rcts_prds, params):
    """ This function finds duplicates within the list of
        reactants and products and converts the params
        to a list of tuples of tuples.

        :param rcts_prds: reaction keys (i.e., reactant/product sets)
        :type rcts_prds: list of tuples
            [((rct1, rct2,...),(prd1, prd2,...),(third body,)),]
        :param params: reaction parameters
        :type params: list [(params1),(params2)]
        :return unique_list: unique list of reactants and products
        :rtype: list of tuples
        [((rct1, rct2,...),(prd1, prd2,...),(third body,)),]
        :return unique_params: updated reaction parameters
            with duplicates included
        :rtype: list of tuples [((params1),(params1')),...]
    """
    # Get the unique entries and number of occurrences
    unique_list = list(set(rcts_prds))
    val = collections.Counter(rcts_prds)
    unique_params = []
    # Loop through each unique item, looking for duplicates
    for idx, rxn in enumerate(unique_list):

        # allocate first value in list
        first = rcts_prds.index(unique_list[idx])
        params_all = [params[first]]

        for i in range(1, val[rxn]):
            first_plus_i = rcts_prds.index(unique_list[idx], first+i)
            params_all.append(params[first_plus_i])

        # check for duplicate plogs
        params_all = split_plog_dct(params_all)
        # convert to tuple and overwrite the parameters
        unique_params.append(tuple(params_all))

    return unique_list, unique_params


def split_plog_dct(param_dct_entry):
    """ Splits 2 or more PLOG dictionaries if they share the same pressures

        :param_dct_entry: values of one rxn_param_dct
        :type tuple(tuple)
        :return param_dct
        :rtype tuple(tuple)
    """
    # convert param_dct [(tup)] to [[lst]]
    param_dct_lst = [list(param_dct_i) for param_dct_i in param_dct_entry]
    # extract plog dictionaries
    plog = np.array(
        [param_dct_vals[4] is not None for param_dct_vals in param_dct_lst],
        dtype=int)
    mask_plog = np.where(plog == 1)[0]

    for plog_i in mask_plog:
        param_dct_vals = param_dct_lst[plog_i]
        keys, plog_params = zip(*list(param_dct_vals[4].items()))

        if any((int(len(param_i)/3) > 1 for param_i in plog_params)):
            new_sets = []

            for idx in range(max([int(len(param_i)/3)
                                  for param_i in plog_params])):
                param_dct_vals_idx = copy.deepcopy(param_dct_vals)
                param_dct_vals_idx[4] = {}
                keys_idx = np.array(
                    [int(len(par_i)/3/(idx+1)) >= 1 for par_i in plog_params],
                    dtype=int)
                mask_keys = np.where(keys_idx == 1)[0]
                new_keys = [keys[i] for i in mask_keys]
                for key_i in new_keys:
                    param_dct_vals_idx[4][key_i] = (
                        param_dct_vals[4][key_i][3 * idx:3*(idx+1)])
                # append the new dictionary to param_dct_lst
                new_sets.append(param_dct_vals_idx)
            # replace the original dictionary
            param_dct_lst[plog_i] = new_sets[0]
            # add new sets
            param_dct_lst.extend(new_sets[1:])

    # convert back to tuples
    param_dct_entry = [tuple(param_dct_i) for param_dct_i in param_dct_lst]
    return param_dct_entry


def get_ea_conv_factor(ea_units):
    """ Get the factor for converting Ea to the desired units of kcal/mole

        :param ea_units: units of activation energies
        :type ea_units: string
        :return ea_conv_factor: conversion factor from ea_units to cal/mol
        :rtype: float
    """
    if ea_units == 'cal/mole':
        ea_conv_factor = 1
    elif ea_units == 'kcal/mole':
        ea_conv_factor = phycon.KCAL2CAL
    elif ea_units == 'joules/mole':
        ea_conv_factor = phycon.J2CAL
    elif ea_units == 'kjoules/mole':
        ea_conv_factor = phycon.KJ2CAL
    elif ea_units == 'kelvins':
        ea_conv_factor = phycon.KEL2CAL
    else:
        raise NotImplementedError(
            f"Invalid ea_units: {ea_units}." +
            " Options: 'kcal/mole', 'cal/mole', 'joules/mole'," +
            " 'kjoules/mole', 'kelvins'"
        )

    return ea_conv_factor


def get_a_conv_factor(rxn_dstr, a_units):
    """ Get the factor for converting A to the desired basis of moles

        :param rxn_dstr: data string for species in reaction block
        :param a_units: units of pre-exponential factor
        :type rxn_dstr, a_units: string
        :return a_conv_factor: conversion factor from a_units to moles
        :rtype: float
    """
    # Get the molecularity
    rcts = reactant_names(rxn_dstr)
    if not isinstance(rcts, tuple):  # convert to list to avoid mistake
        rcts = [rcts]
    molecularity = len(rcts)

    # Find out whether there is a third body
    trd_body = third_body(rxn_dstr)[0]
    # if 3rd body has '(', no effect on units
    if trd_body is not None and '(+' not in trd_body:
        molecularity += 1

    if a_units == 'moles':
        a_conv_factor = 1
    elif a_units == 'molecules':
        a_conv_factor = phycon.NAVO**(molecularity-1)
    else:
        raise NotImplementedError(
            f"Invalid a_units: {a_units}. Options: 'moles' or 'molecules'"
        )

    return a_conv_factor


# Not used here, but called elsewhere
def ratek_fit_info(rxn_dstr):
    """ Read the information describing features of the fits to the
        rate constants
    """

    # Read the temperatures and the Errors from the lines
    pressure_ptt = (
        'Pressure:' + app.SPACES +
        app.capturing(app.one_of_these([app.FLOAT, 'High']))
    )
    trange_ptt = (
        'Temps: ' + app.SPACES +
        app.capturing(app.INTEGER) + '-' + app.capturing(app.INTEGER) +
        app.SPACES + 'K'
    )
    mean_ptt = (
        'MeanAbsErr:' + app.SPACES +
        app.capturing(app.FLOAT) + app.escape('%') +
        ','
    )
    max_ptt = (
        'MaxErr:' + app.SPACES +
        app.capturing(app.FLOAT) + app.escape('%')
    )
    pressure_caps = apf.all_captures(pressure_ptt, rxn_dstr)
    trange_caps = apf.all_captures(trange_ptt, rxn_dstr)
    mean_caps = apf.all_captures(mean_ptt, rxn_dstr)
    max_caps = apf.all_captures(max_ptt, rxn_dstr)

    pressures = []
    for pressure in pressure_caps:
        if pressure != 'High':
            pressures.append(float(pressure))
        elif pressure == 'High':
            pressures.append(pressure)
    trange_vals = []
    for cap in trange_caps:
        temp1, temp2 = cap
        trange_vals.append([int(temp1), int(temp2)])
    if mean_caps is not None:
        mean_vals = [float(val) for val in mean_caps]
    else:
        mean_vals = []
    if max_caps is not None:
        max_vals = [float(val) for val in max_caps]
    else:
        max_vals = []

    # Build the inf_dct
    inf_dct = {}
    for idx, pressure in enumerate(pressures):
        inf_dct[pressure] = {'temps': trange_vals[idx]}
        if mean_vals:
            inf_dct[pressure].update({'mean_err': mean_vals[idx]})
        if max_vals:
            inf_dct[pressure].update({'max_err': max_vals[idx]})

    return inf_dct
