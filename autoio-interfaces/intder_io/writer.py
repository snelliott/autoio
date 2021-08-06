""" Write the INTDER main and auxiliary input file
"""

import intder_io._format as intder_format


def input_file(geo, zma):
    """ Write the main INTDER input file. Currently
        just supports a basic harmonic frequency and
        total energy distribution calculation.

        zma and geo must align

        :param geo: geometry to build input for
        :type geo: automol geometry data structure
        :param zma: Z-Matrix corresponding to geometry
        :type zma: automol Z-Matrix data structure
        :rtype: str
    """

    inp_str = (
        intder_format.header_format(geo) + '\n' +
        intder_format.internals_format(zma) + '\n' +
        intder_format.geometry_format(geo) + '\n' +
        intder_format.symbols_format(geo)
    )

    return inp_str


def cart_hess_file(hess):
    """ Write a file with the Cartesian Hessian auxiliary
        input that corresponds to the FCMINT file for CFOUR.

        :param hess: mass-weighted Hessian (in a.u.)
        :type hess: numpy.ndarray
        :rtype: str
    """

    natom = len(hess) // 3
    hess_str = '{0:>5d}{1:>5d}\n'.format(natom, natom*3)
    hess_str += intder_format.hessian_format(hess)

    return hess_str
