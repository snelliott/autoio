""" Functions that take data parsed from INTDER output
    as well as other functionality and calculate useful
    values
"""

from automol.zmat import distance_coordinate_name
from automol.zmat import central_angle_coordinate_name
from automol.zmat import dihedral_angle_coordinate_name


def ted_zmatrix_coordinates(zma, mode_idx, intl_coords, ted_assignments):
    """ Take the output of a Total Energy Distribution analysis
        from INTDER and determine what internal coordinates
        defined in a Z-Matrix make up that mode.

        The output has the internal coordinates defined in terms of
        atom indices of a Cartesian geometry+Hessian. Therefore, the
        input zmatrix must have an atom ordering that aligns with
        that Cartesian system used to generate the TED analysis.

        See intder_io.reader documentation for what intl_coords and
        ted_assignment variable objects are structured like.

        :param zma: Z-Matrix with coordinates to get names for
        :typ zma: automol.zmat object
        :param intl_coords: definition of internal coordinates in INTDER
            in terms of the atom indices
        :type inl_coords: tuple(tuple(str, tuple(int)))
    """

    # Get the intl coord idxs associated with vibrational mode
    ted_mode = ted_assignments[mode_idx]
    ted_mode_intls = tuple(ted_mode[1].keys())

    # For each internal coordinate for the TED,
    # Get the type and atom and use this to get zmat coordinate name
    ted_zmat_names = ()
    for intl_idx in ted_mode_intls:
        typ, idxs = intl_coords[abs(intl_idx)]
        print(typ, idxs)
        if typ == 'STRE':
            name = distance_coordinate_name(zma, *idxs)
        elif typ == 'BEND':
            name = central_angle_coordinate_name(zma, *idxs)
        elif typ == 'TORS':
            name = dihedral_angle_coordinate_name(zma, *idxs)
        print(name)
        ted_zmat_names += (name,)

    return ted_zmat_names
