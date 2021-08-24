""" write the code

    need a way to write the fitting errors
"""

import chemkin_io
import autoreact


def chemkin_mechanism(rxn_param_dct):
    """ Writes the chemkin strings that contain the
        fitting parameters for all reactions in the
        parameter dictionary.
    """

    ckin_str = ''
    for rxn, param_obj in rxn_param_dct.items():

        # Generate the string for the reactions and products
        rxn_str = chemkin_io.writer.format_rxn_name(rxn)

        # Write the parameters to the string
        param_str = chemkin_reaction(rxn_str, param_obj, max_length=45)

        # Write the fitting errors (NEED object to do this)
        fit_str = '! FIT ERR STRING'

        # Write to overall string
        ckin_str += param_str + fit_str + '\n\n'

    return ckin_str


def chemkin_reaction(reaction, params, max_length=45):
    """ Write a chemkin string from a params object.

        :param reaction: reaction
        :type reaction: str
        :param params: parameters for a reaction
        :type params: autoreact.ReactionParams object
        :param max_length: length of the longest reaction name in the mechanism
        :type max_length: int
        :param mex
        :rtype: str
    """

    # Determine which functional form to write to string
    # Based on existing parameters and hierarchy of expressions
    form = params.form_to_write()

    # Write the string for the selected functional form
    if form == 'arrhenius':
        ckin_str = chemkin_io.writer.reaction.arrhenius(
            reaction,
            params.arrhenius,
            max_length=max_length)
    elif form == 'plog':
        ckin_str = chemkin_io.writer.reaction.plog(
            reaction,
            params.plog,
            max_length=max_length)
    elif form == 'chebyshev':
        ckin_str = chemkin_io.writer.reaction.chebyshev(
            reaction,
            params.chebyshev['highp_arr'],
            params.chebyshev['alpha'],
            params.chebyshev['tlim'][0],
            params.chebyshev['tlim'][1],
            params.chebyshev['plim'][0],
            params.chebyshev['plim'][1],
            max_length=max_length)
    elif form == 'troe':
        ckin_str = chemkin_io.writer.reaction.troe(
            reaction,
            params.troe['highp_arr'],
            params.troe['lowp_arr'],
            params.troe['troe_params'],
            colliders=params.colliders,
            max_length=max_length)
    elif form == 'lindemann':
        ckin_str = chemkin_io.writer.reaction.lindemann(
            reaction,
            params.lindemann['highp_arr'],
            params.lindemann['lowp_arr'],
            colliders=params.colliders,
            max_length=max_length)

    else:
        ckin_str = None

    return ckin_str


if __name__ == '__main__':

    # Make reactions
    _r1 = (('W1',), ('P1',), (None,))
    _p1 = autoreact.params.RxnParams()
    _p1.set_arr([[1e12, 1.5, 50000], [2e12, 2.5, 50000]])

    _r2 = (('W1',), ('P2',), (None,))
    _p2 = autoreact.params.RxnParams()
    _p2.set_plog({'high': [[1e14, 2.5, 60000]],
                  1.0: [[2e14, 2.5, 60000]],
                  10.0: [[2e14, 2.5, 60000]]})

    _r3 = (('W1',), ('P3',), (None,))
    _p3 = autoreact.params.RxnParams()
    _p3.set_chebyshev([[1e16, 3.5, 70000]], [[1e12, 1.8, 45000]],
                      [100.0, 500.0], [1.0, 100.0],
                      [[1.0, 2.0, 3.0, 4.0], [1.1, 2.1, 3.1, 4.1],
                       [1.3, 2.3, 3.3, 4.3], [1.4, 2.4, 3.4, 4.4]])

    _r4 = (('W1',), ('P4',), (None,))
    _p4 = autoreact.params.RxnParams()
    _p4.set_troe([[1e16, 3.5, 70000]], [[1e12, 1.8, 45000]],
                 55.0, [1.0, 2.0, 3.0, 4.0])

    _r5 = (('W1',), ('P5',), (None,))
    _p5 = autoreact.params.RxnParams()
    _p5.set_lindemann([[1e16, 3.5, 70000]], [[1e12, 1.8, 45000]])

    # Make a parameter dictionary
    _rxn_param_dct = {
        _r1: _p1,
        _r2: _p2,
        _r3: _p3,
        _r4: _p4,
        _r5: _p5
    }

    # Run code
    _ckin_str = chemkin_mechanism(_rxn_param_dct)
    print(_ckin_str)
