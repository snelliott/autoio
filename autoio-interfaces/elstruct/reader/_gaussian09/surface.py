""" potential energy surface information readers
"""

import numpy
from phydat import ptab
import automol
import autoread as ar
from autoparse import cast as _cast
import autoparse.pattern as app
import autoparse.find as apf


def gradient(output_str):
    """ Reads the molecular gradient (in Cartesian coordinates) from
        the output file string. Returns the gradient in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    grad = ar.matrix.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape('Forces (Hartrees/Bohr)'), app.NONNEWLINE),
            app.LINE, app.LINE, '']),
        line_start_ptt=app.LINESPACES.join([app.UNSIGNED_INTEGER] * 2))
    if grad is not None:
        grad = numpy.multiply(grad, -1.0)

    return grad


def hessian(output_str):
    """ Reads the molecular Hessian (in Cartesian coordinates) from
        the output file string. Returns the Hessian in atomic units.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    comp_ptt = app.UNSIGNED_INTEGER
    mat = ar.matrix.read(
        output_str,
        val_ptt=app.EXPONENTIAL_FLOAT_D,
        start_ptt=(
            app.escape('Force constants in Cartesian coordinates:') +
            app.lpadded(app.NEWLINE)),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=True)

    if mat is not None:
        mat = [[_cast(apf.replace('d', 'e', dst, case=False)) for dst in row]
               for row in mat]

    if mat is None:
        comp_ptt = app.one_of_these(['X', 'Y', 'Z']) + app.UNSIGNED_INTEGER
        mat = ar.matrix.read(
            output_str,
            start_ptt=(app.escape('The second derivative matrix:') +
                       app.lpadded(app.NEWLINE)),
            block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                             app.padded(app.NEWLINE)),
            line_start_ptt=comp_ptt,
            tril=True)

    if mat is not None:
        mat = tuple(map(tuple, mat))

    return mat


def harmonic_frequencies(output_str):
    """ Reads the harmonic vibrational frequencies from
        the output file string. Returns the frequencies in cm-1.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(float)
    """

    pattern = 'Frequencies --' + app.capturing(app.LINE_FILL)
    captures = apf.all_captures(pattern, output_str)
    if captures is not None:
        freqs = []
        for capture in captures:
            vals = capture.split()
            for val in vals:
                freqs.append(float(val))
    else:
        freqs = None
    return freqs


def normal_coordinates(output_str):
    """ Reads the displacement along the normal modes (in Cartesian coordinates)
        from the output string. Returns the coordinates in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(tuple(float))
    """

    # Set the patterns to read the normal coordinates
    comp_ptt = app.UNSIGNED_INTEGER + app.SPACES + app.UNSIGNED_INTEGER

    # Gaussian prints freqs/norm coords as columns. Up to 3 columns printed
    # Three patterns used to handle cases of when 1, 2, or 3 columns printed
    start_ptt = 'Atom  AN      X      Y      Z'
    start_ptt2 = '        X      Y      Z'
    start_ptt_lst = (
        (start_ptt + start_ptt2 + start_ptt2),  # 3 freq column(s)
        (start_ptt + start_ptt2),               # 2 freq column(s)
        (start_ptt),                            # 1 freq column(s)
    )

    # Read the normal coordinate modes
    nmats = []
    for mode in apf.split('Frequencies', output_str)[1:]:
        mat = None
        for start in start_ptt_lst:
            _ptt = app.padded(app.NEWLINE).join([app.escape(start), ''])
            mat = ar.matrix.read(
                mode,
                start_ptt=_ptt,
                line_start_ptt=comp_ptt)
            if mat is not None:
                break
        nmats.append(mat)

    # Parse the coordinates out of the matrices if all matrices found
    if all(mat is not None for mat in nmats):
        nmodes = ()
        for mat in nmats:
            nmat = numpy.array(mat)
            for i in range(int(len(nmat)/3)):
                nmodes += (nmat[:, i*3:(i+1)*3],)
    else:
        nmodes = None

    print('normal coords test')
    print(nmodes)

    return nmodes


def irc_points(output_str):
    """ Reads the geometries, gradients, and Hessians at each point along the
        Intrinsic Reaction Coordinate from the output string.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: (geom data structure, tuple(tuple(float)), tuple(tuple(float)))
    """

    # Lines
    output_lines = output_str.splitlines()

    # Obtain all of the info for the points off the saddle point
    # Find the lines with point number to get the strings
    section_starts = []
    for i, line in enumerate(output_lines):
        if 'Point Number  1 in' in line:
            section_starts.append(i)
            break
    for i, line in enumerate(output_lines):
        if 'OF POINTS ALONG THE PATH' in line:
            section_starts.append(i)

    # Get the TS info first
    sadpt_str = '\n'.join(output_lines[0:section_starts[0]])
    sadpt_geom = sadpt_geometry(sadpt_str)
    # sadpt_grad = gradient(sadpt_str)
    # sadpt_hess = hessian(sadpt_str)

    # Now start getting other points by getting a list of each string with info
    pt_strs = []
    for i in range(1, len(section_starts)):
        start_line = section_starts[i-1]
        end_line = section_starts[i]
        pt_strs.append('\n'.join(output_lines[start_line+1:end_line]))

    # Obtain the grads and hessians
    geoms = []
    grads = []
    hessians = []
    for string in pt_strs:
        geoms.append(irc_geometry(string))
        # pt_grad = gradient(string)
        # pt_hess = hessian(string)
        # if pt_grad is not None:
        #     grads.append(pt_grad)
        # if pt_hess is not None:
        #     hessians.append(pt_hess)

    # Combine with the 0 index info
    geoms = [sadpt_geom] + geoms
    if grads:
        grads = [] + grads
        # grads = [sadpt_grad] + grads
    if hessians:
        hessians = [] + hessians
        # hessians = [sadpt_hess] + hessians

    return geoms, grads, hessians


def sadpt_geometry(sadpt_str):
    """ Reads the molecular geometry of the saddle point (first point on the
        Intrinsic Reaction Coordinate). Returns the geometry in Bohr.

        :param sadpt_str: string of the output file containing the saddle point
        :type sadpt_str: str
        :rtype: automol geom data structure
    """

    nums, xyzs = ar.geom.read(
        sadpt_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('Input orientation:'),
            app.LINE, app.LINE, app.LINE, app.LINE, '']),
        symb_ptt=app.UNSIGNED_INTEGER,
        line_start_ptt=app.UNSIGNED_INTEGER,
        line_sep_ptt=app.UNSIGNED_INTEGER)
    syms = tuple(map(ptab.to_symbol, nums))
    geo = automol.geom.from_data(syms, xyzs, angstrom=True)

    return geo


def irc_geometry(output_str):
    """ Reads the molecular geometry at a point on the
        Intrinsic Reaction Coordinate. Returns the geometry in Bohr.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(automol geom data structure)
    """

    nums, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('CURRENT STRUCTURE'),
            app.LINE, app.LINE, app.LINE, app.LINE, app.LINE, '']),
        symb_ptt=app.UNSIGNED_INTEGER,
        line_start_ptt=app.UNSIGNED_INTEGER)
    if any(val is None for val in (nums, xyzs)):
        nums, xyzs = ar.geom.read(
            output_str,
            start_ptt=app.padded(app.NEWLINE).join([
                app.escape('Input orientation:'),
                app.LINE, app.LINE, app.LINE, app.LINE, '']),
            symb_ptt=app.UNSIGNED_INTEGER,
            line_start_ptt=app.UNSIGNED_INTEGER,
            line_sep_ptt=app.UNSIGNED_INTEGER,)

    if all(val is not None for val in (nums, xyzs)):
        syms = tuple(map(ptab.to_symbol, nums))
        geo = automol.geom.from_data(syms, xyzs, angstrom=True)
    else:
        geo = None

    return geo


def irc_path(output_str):
    """ Reads the coordinates and electronic energies (relative to saddple point)
        of the Intrinsic Reaction Coordinate summarized at the end of
        the output file.
        Returns the energy in Hartress.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: tuple(automol geom data structure)
    """

    # Reads the coordiantes
    coordinates = _read_irc_reaction_path_summary(output_str, 'coord')

    # Reads the energies (the ts/sadpt)
    ptt = (
        'Energies reported relative to the TS energy of' +
        app.SPACES +
        app.capturing(app.FLOAT)
    )
    ts_energy = apf.last_capture(ptt, output_str)
    pt_energies = _read_irc_reaction_path_summary(output_str, 'energy')
    if ts_energy and pt_energies:
        energies = [float(ts_energy) + ene for ene in pt_energies]

    # See if the enes need to be flipped so the ts ene is first
    if pt_energies[0] != 0.0:
        coordinates = coordinates[::-1]
        energies = energies[::-1]

    return (coordinates, energies)


def _read_irc_reaction_path_summary(output_str, read_val):
    """ Reads the values for the Intrinsic Reaction Path from the table.

        :param output_str: string of the program's output file
        :type output_str: str
        :param read_val: value to read from table
        :type read_val: str
        :rtype: tuple(automol geom data structure)
    """

    assert read_val in ('energy', 'coord')

    block = apf.last_capture(
        (app.escape('Summary of reaction path following') +
         app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
         app.escape('Total number of points:') + app.SPACES + app.INTEGER),
        output_str)

    if read_val == 'energy':
        pattern = (
            app.INTEGER + app.SPACES +
            app.capturing(app.FLOAT) +
            app.SPACES +
            app.FLOAT
        )
    elif read_val == 'coord':
        pattern = (
            app.INTEGER + app.SPACES +
            app.FLOAT +
            app.SPACES +
            app.capturing(app.FLOAT)
        )

    captures = apf.all_captures(pattern, block)
    if captures is not None:
        values = [float(capture) for capture in captures]
    else:
        values = None

    return values


if __name__ == '__main__':
    with open('output.dat', 'r') as fobj:
        _str = fobj.read()
    normal_coordinates(_str)
