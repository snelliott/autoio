""" Functions to format the data into strings appropriate
"""

import os
from itertools import chain
import automol.geom
import automol.zmat
import automol.util
from ioformat import build_mako_str
from ioformat import indent


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def header_format(geo):
    """ format the header
    """

    natom = automol.geom.count(geo)
    if not automol.geom.is_linear(geo):
        nintl = 3 * natom - 6
    else:
        nintl = 3 * natom - 5

    header_keys = {
        'natom': natom,
        'nintl': nintl,
        'comment': 'TED Calculation'
    }

    return build_mako_str(
        template_file_name='header.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=header_keys,
        remove_whitespace=False)


def internals_format(zma):
    """ Format the strings with the internal coordinates
    """

    key_mat = automol.zmat.key_matrix(zma)

    # Write the stretch, bend, and torsion coordinates
    intl_str = ''
    for i, row in enumerate(key_mat[1:]):
        intl_str += 'STRE{0:>6d}{1:>5d}\n'.format(
            i+1+1, row[0]+1)
    for i, row in enumerate(key_mat[2:]):
        intl_str += 'BEND{0:>6d}{1:>5d}{2:>5d}\n'.format(
            i+2+1, row[0]+1, row[1]+1)
    for i, row in enumerate(key_mat[3:]):
        intl_str += 'TORS{0:>6d}{1:>5d}{2:>5d}{3:>5d}\n'.format(
            i+3+1, row[0]+1, row[1]+1, row[2]+1)

    # Remove final newline character
    intl_str = intl_str.rstrip()

    return intl_str


def geometry_format(geo):
    """ Format the geometry string
    """

    # Build geom str
    geo_str = ''
    for (_, xyz) in geo:
        geo_str += '{:>14.5f}{:>14.5f}{:>14.5f}\n'.format(*xyz)

    # Indent the lines and remove final newline character
    geo_str = indent(geo_str, 4)
    geo_str = geo_str.rstrip()

    return geo_str


def symbols_format(geo):
    """ Format the symbols
    """

    symbs = automol.geom.symbols(geo)
    symb_str = automol.util.vec.string(
        symbs, num_per_row=6, val_format='{0:>12s}')

    return symb_str


def hessian_format(hess):
    """ Format a mass-weighted Hessian into a string for the
        auxiliary input file for INTDER.

        :param hess: mass-weighted Hessian
        :type hess: numpy.ndarray
        :rtype: str
    """

    # Flatten the Hessian out for ease of writing the string
    hess = tuple(chain.from_iterable(hess))
    hess_str = automol.util.vec.string(
        hess, num_per_row=3, val_format='{0:>20.10f}')

    return hess_str
