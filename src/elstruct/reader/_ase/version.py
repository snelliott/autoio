""" Reads the program name and version number from the output file
"""

import autoparse.pattern as app
import autoparse.find as apf


def program_name(output_dct):
    """ Reads the program name from the output dictionary.
        :param output_dct: a dictionary of output information
        :type output_dct: dict
        :rtype: str
    """
    
    return output_dct['version'].split('_')[0].lower()


def program_version(output_dct):
    """ Reads the program version number from the output dictionary.

        :param output_dct: a dictionary of output information
        :type output_dct: dict
        :rtype: str
    """

    return output_dct['version']
