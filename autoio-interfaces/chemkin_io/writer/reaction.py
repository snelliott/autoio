"""
Writes Chemkin-formatted strings containing the rate parameters
"""

from chemkin_io.writer import _util

# Define the name_buffer between the longest reaction name and the Arr params
BUFFER = 5


def get_ckin_str(rxn, params, max_len=45):
    """ Writes a reaction string from a RxnParams object

        :param rxn: tuple describing reactants, products, and third body
        :type rxn: tuple ((rct1, rct2, ...), (prd1, prd2, ...), (thirdbod,))
        :param params: parameters for a reaction
        :type params: autoreact.RxnParams object
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :return ckin_str: Chemkin-formatted string describing the reaction
        :rtype: str
    """

    # Get the Chemkin string for the chemical reaction
    reaction = _util.format_rxn_name(rxn)

    # Determine which functional form to write to string
    form = params.form_to_write()

    # Write the string for the selected functional form
    if form == 'arr':
        ckin_str = arr(reaction, params.arr, colliders=params.arr_collid, 
                       max_len=max_len)
    elif form == 'plog':
        ckin_str = plog(reaction, params.plog, max_len=max_len)
    elif form == 'cheb':
        ckin_str = cheb(reaction, params.cheb['alpha'], params.cheb['tlim'], 
                        params.cheb['plim'], params.cheb['one_atm_arr'],
                        max_len=max_len)
    elif form == 'troe':
        ckin_str = chemkin_io.writer.reaction.troe(
            reaction,
            params.troe['highp_arr'],
            params.troe['lowp_arr'],
            params.troe['troe_params'],
            colliders=params.colliders,
            max_len=max_len)
    elif form == 'lindemann':
        ckin_str = chemkin_io.writer.reaction.lindemann(
            reaction,
            params.lindemann['highp_arr'],
            params.lindemann['lowp_arr'],
            colliders=params.colliders,
            max_len=max_len)

    else:
        ckin_str = None

    return ckin_str


def troe(reaction, high_params, low_params, troe_params, colliders=None, max_len=45):
    """ Writes a reaction in the Troe form

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param high_params: Arrhenius high-P parameters
        :type high_params: list of floats
        :param low_params: Arrhenius low-P parameters
        :type low_params: list of floats
        :param troe_params: Troe parameters: alpha, T***, T*, and T**
            (T** is optional)
        :type troe_params: list of floats
        :param colliders: collision enhancement factors for bath gases
        :type colliders: list((str, float))
        :return troe_str: Chemkin reaction string with Troe parameters
        :rtype: str
    """
    assert len(high_params[0]) == 3, (
        f'{len(high_params)} highP params for {reaction}, should be 3'
    )
    assert len(low_params[0]) == 3, (
        f'{len(low_params)} lowP params for {reaction}, should be 3'
    )
    assert len(troe_params) in (3, 4), (
        f'{len(troe_params)} Troe params for {reaction}, should be 3 or 4'
    )

    troe_str = _highp_str(reaction, high_params[0], max_len=max_len)
    troe_str += _lowp_str(low_params[0], max_len=max_len)
    troe_str += _misc_troe_and_cheb('TROE', troe_params, newline=True, 
                                    val='exp')

    # Write the collider efficiencies string
    if colliders:
        troe_str += _format_collider_string(colliders)

    return troe_str


def lind(reaction, high_params, low_params, colliders=None, max_len=45):
    """ Writes a reaction in the Lindemann form

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param high_params: Arrhenius high-P parameters
        :type high_params: list of floats
        :param low_params: Arrhenius low-P parameters
        :type low_params: list of floats
        :param colliders: names and collision enhancement factors for bathgases
        :type colliders: list((str, float))
        :return lind_str: Chemkin reaction string with Lindemann parameters
        :rtype: str
    """
    lind_str = _highp_str(reaction, high_params[0], max_len=max_len)
    lind_str += _lowp_str(low_params[0], max_len=max_len)

    # Write the collider efficiencies string
    if colliders:
        lind_str += _format_collider_string(colliders)

    return lind_str


def plog(reaction, plog_param_dct, max_len=45):
    """ Writes a reaction in the PLOG form

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param plog_param_dct: Arrhenius fitting parameters at all pressures
        :type plog_param_dct: dict{pressure: [Arrhenius params]}
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :return plog_str: Chemkin reaction string with PLOG parameters
        :rtype: str
    """

    def _pressure_str(pressure, params, max_len=45):
        """ Write a line in a PLOG string

            :param pressure: pressure at which to write the string
            :type pressure: float
            :param params: Arrhenius parameters at the specified pressure
            :type params: list of floats
            :param max_len: length of the longest reaction name
                in the mechanism
            :type max_len: int
            :return single_str: Chemkin reaction string with PLOG parameters
                at a single pressure
            :rtype: str
        """

        plog_buffer = str(max_len + BUFFER - 12)
        [a_par, n_par, ea_par] = params
        single_str = (
            '{0:<' + plog_buffer +
            's}{1:<12.3E}{2:<10.3E}{3:>9.3f}{4:>9.0f} /\n').format(
                '    PLOG /', pressure, a_par, n_par, ea_par)

        return single_str

    # Obtain a list of the pressures and sort from low to high pressure
    unsorted_pressures = plog_param_dct.keys()
    pressures = sorted(unsorted_pressures)

    # Write the header for the reaction, which includes the 1-atm fit if avail
    if 1 in pressures:
        if len(plog_param_dct[1]) > 1:
            comment = (
                'Duplicates exist at 1 atm (see below);' +
                ' only single 1-atm fit is written')
            plog_str = _highp_str(
                reaction, plog_param_dct[1][0], max_len=max_len,
                inline_comment=comment)
        else:
            comment = 'Arrhenius parameters at 1 atm'
            plog_str = _highp_str(
                reaction, plog_param_dct[1][0], max_len=max_len,
                inline_comment=comment)
    else:
        comment = 'No fit at 1 atm available'
        plog_str = _highp_str(reaction, [1.0, 0.0, 0.0], max_len=max_len,
                              inline_comment=comment)

    # Loop over each pressure
    for pressure in pressures:
        plog_params = plog_param_dct[pressure]
        for param_set in plog_params:
            assert len(param_set) % 3 == 0, (
                f'Arr params should be a multiple of 3, is {len(param_set)}' +
                f' for {reaction}'
            )

            # Loop over however many Arrhenius sets there are,
            # writing a PLOG line for each
            num_arr_sets = int(len(param_set)/3)
            for idx in range(num_arr_sets):
                current_param_set = param_set[3*idx:3*(idx+1)]
                plog_str += _pressure_str(
                    pressure, current_param_set, max_len=max_len)

    return plog_str


def cheb(reaction, alpha, tlim, plim, one_atm_arr=None, max_len=45):
    """ Writes a reaction in the Chebyshev form

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param alpha: Chebyshev coefficient matrix
        :type alpha: numpy.ndarray
        :param tlim: Chebyshev temperature limits
        :type tlim: tuple (tmin, tmax)
        :param plim: Chebyshev pressure limits
        :type plim: tuple (pmin, pmax)
        :param one_atm_arr: Arrhenius parameters at 1 atm
        :type one_atm_arr: tuple of tuples ((A1, n1, Ea1), (A2, n2, Ea2), ...)
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :return cheb_str: Chemkin reaction string with Chebyshev parameters
        :rtype: str
    """

    # Write reaction header (with third body added) and high-pressure params
    if one_atm_arr is not None:
        comment = 'Arrhenius parameters at 1 atm'
    else:  # if the params look fake
        comment = None
        one_atm_arr = [[1, 0, 0],]
    cheb_str = _highp_str(reaction, one_atm_arr[0], max_len=max_len,
                          inline_comment=comment)

    # Write the temperature and pressure ranges
    cheb_str += _misc_troe_and_cheb(
        'TCHEB', tlim, newline=True, val='float')
    cheb_str += _misc_troe_and_cheb(
        'PCHEB', plim, newline=True, val='float')

    # Write the dimensions of the alpha matrix
    nrows = len(alpha)
    ncols = len(alpha[0])
    cheb_str += _misc_troe_and_cheb(
        'CHEB', (nrows, ncols), newline=True, val='int')

    # Write the parameters from the alpha matrix
    for row in alpha:
        cheb_str += _misc_troe_and_cheb('CHEB', row, newline=True, val='exp')

    return cheb_str


def arr(reaction, arr_tuples, colliders=None, max_len=45):
    """ Writes a reaction in the Arrhenius form

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param arr_tuples: Arrhenius high-P (i.e., high-P) parameters
        :type arr_tuples: tuple of tuples ((A1, n1, Ea1), (A2, n2, Ea2), ...)
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :return arr_str: Chemkin reaction string with Arrhenius parameters
        :rtype: str
    """

    # Write the Arrhenius parameter string
    arr_str = ''
    for arr_tuple in arr_tuples:
        arr_str += _highp_str(reaction, arr_tuple, max_len=max_len)
        if len(arr_tuples) > 1:
            arr_str += 'DUP\n'

    # Write the collider efficiencies string
    if colliders is not None:
        arr_str += _format_collider_string(colliders)

    return arr_str


# Various formatting functions
def fit_info(pressures, temp_dct, err_dct):
    """ Write the string detailing the temperature ranges and fitting errors
        associated with the rate-constant fits at each pressure.

        :param pressures: pressures the k(T,P)s were calculated at
        :type pressures: list(float)
        :param temp_dct: temperature ranges (K) fits were done at each pressure
        :type temp_dct: dict[pressure, [temp range]]
        :param err_dct: errors associated with the fits at each pressure
        :type err_dct: dict[pressure, [mean err, max err]]
        :return inf_str: string containing all of the fitting info
        :rtype: str
    """
    # Make temp, err dcts empty if fxn receives None; add 'high' to pressures
    temp_dct = temp_dct if temp_dct else {}
    err_dct = err_dct if err_dct else {}
    if 'high' in temp_dct or 'high' in err_dct:
        pressures += ['high']

    # Check the temp and err dcts have same presures as rate_dcts
    if temp_dct:
        assert set(pressures) == set(temp_dct.keys())
    err_dct = err_dct if err_dct else {}
    if err_dct:
        assert set(pressures) == set(err_dct.keys())

    # Write string showing the temp fit range and fit errors
    inf_str = '! Info Regarding Rate Constant Fits\n'
    for pressure in pressures:
        if temp_dct:
            [min_temp, max_temp] = temp_dct[pressure]
            temps_str = '{0:.0f}-{1:.0f} K'.format(
                min_temp, max_temp)
            temp_range_str = 'Temps: {0:>12s}, '.format(
                temps_str)
        else:
            temp_range_str = ''
        if err_dct:
            [mean_err, max_err] = err_dct[pressure]
            err_str = '{0:11s} {1:>5.1f}%,  {2:7s} {3:>5.1f}%'.format(
                'MeanAbsErr:', mean_err, 'MaxErr:', max_err)
        else:
            err_str = ''

        # Put together the whole info string
        if pressure != 'high':
            pstr = '{0:<9.3f}'.format(pressure)
        else:
            pstr = '{0:<9s}'.format('High')
        inf_str += '! Pressure: {0} {1} {2}\n'.format(
            pstr, temp_range_str, err_str)

    return inf_str


def _highp_str(reaction, arr_tuple, max_len=45, inline_comment=None):
    """ Write a single set of high-pressure Arrhenius parameters

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param arr_tuple: Arrhenius high-P parameters
        :type arr_tuple: tuple (or list) of floats (A, n, Ea)
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :param inline_comment: optional inline comment
        :type inline_comment: str
        :return highp_str: Chemkin reaction string with high-P Arrhenius params
        :rtype: str
    """

    highp_buffer = str(max_len + BUFFER)
    [a_par, n_par, ea_par] = arr_tuple
    highp_str = (
        '{0:<' + highp_buffer + 's}{1:<10.3E}{2:>9.3f}{3:>9.0f}').format(
            reaction, a_par, n_par, ea_par)

    if inline_comment:
        highp_str += '   ! ' + inline_comment

    highp_str += '\n'

    return highp_str


def _lowp_str(arr_tuple, max_len=45, inline_comment=None):
    """ Write a single set of low-pressure Arrhenius parameters (i.e., w/'LOW')

        :param reaction: Chemkin-formatted reaction name
        :type reaction: str 
        :param arr_tuple: Arrhenius low-P parameters
        :type arr_tuple: tuple (or list) of floats (A, n, Ea)
        :param max_len: length of the longest reaction name in the mechanism
        :type max_len: int
        :param inline_comment: optional inline comment
        :type inline_comment: str
        :return lowp_str: Chemkin reaction string with low-P Arrhenius params
        :rtype: str
    """

    lowp_buffer = str(max_len + BUFFER)
    [a_par, n_par, ea_par] = arr_tuple
    lowp_str = (
        '{0:<' + lowp_buffer + 's}{1:<10.3E}{2:>9.3f}{3:>9.0f}  /').format(
            '    LOW  /', a_par, n_par, ea_par)

    if inline_comment:
        lowp_str += '   ! ' + inline_comment

    lowp_str += '\n'

    return lowp_str


def _misc_troe_and_cheb(header, params, newline=False, val='exp'):
    """ Write a string containing slash-delimited params used for
        Troe and Chebyshev functional forms

        :param header: the header str to be written before the slashes
        :type header: str
        :param params: params to be written
        :type params: list or numpy array or tuple of floats
        :param newline: whether or not to add a new line after the params
        :type newline: Bool
        :param val: indicates how the parameters are to be formatted;
            should be 'exp' or 'float' or 'int'
        :type val: str
        :return params_str: a string containing the formatted params
        :rtype: str
    """

    if val == 'exp':
        val_str = '{0:12.3E}'
    elif val == 'float':
        val_str = '{0:12.2f}'
    elif val == 'int':
        val_str = '{0:12.0f}'

    params_str = '    {0:5s}/'.format(header.upper())
    for param in params:
        if param:  # skip any None entries
            params_str += ''.join(val_str.format(param))
    params_str += ' /'
    if newline:
        params_str += '\n'

    return params_str


def _format_collider_string(colliders):
    """ Write the string for the bath gas collider and their efficiencies
        for the Lindemann and Troe functional expressions:

        :param colliders: the {collider: efficiency} dct
        :type colliders: dct {str: float}
        :return: collider_str: Chemkin string with colliders and efficiencies
        :rtype: str
    """

    collider_str = ' ' * BUFFER
    for collider, efficiency in colliders.items():
        collider_str += ''.join( ('{0:s}/{1:4.3f}/   '.format(
            collider, efficiency)))
    collider_str += '\n'

    return collider_str
