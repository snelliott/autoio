""" status checkers for ASE output
"""

import elstruct.par

def has_normal_exit_message(output_dct):
    """ Assess whether the output dictionary demonstrates
        normal program termination.

        :param output_dct: a dictionary of output information
        :type output_dct: dict
        :rtype: bool
    """
    return output_dct.get('normal_termination', False)


# Parsers for convergence success messages
def _has_scf_convergence_message(output_dct):
    """ Assess whether the output dictionary indicates
        successful convergence of the SCF procedure.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE calculations will have already failed if SCF didn't converge
    # so if we have any output, SCF converged
    return True


def _has_opt_convergence_message(output_dct):
    """ Assess whether the output dictionary indicates
        successful convergence of the geometry optimization.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE optimizations will have failed if they didn't converge
    # so if we have any output, optimization converged
    return True


def _has_irc_convergence_message(output_dct):
    """ Assess whether the output dictionary indicates
        successful convergence of the IRC calculation.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE IRC calculations will fail if they don't converge
    # so if we have any output, IRC converged
    return True


# Parsers for various error messages
def _has_scf_nonconvergence_error_message(output_dct):
    """ Assess whether the output dictionary indicates
        unsuccessful convergence of the SCF procedure.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE will raise an error if SCF doesn't converge, so we won't get here
    return False


def _has_opt_nonconvergence_error_message(output_dct):
    """ Assess whether the output dictionary indicates
        unsuccessful convergence of the geometry optimization.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE will raise an error if optimization doesn't converge
    return False


def _has_irc_nonconvergence_error_message(output_dct):
    """ Assess whether the output dictionary indicates
        unsuccessful convergence of the IRC calculation.

        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    # ASE will raise an error if IRC doesn't converge
    return False


def _has_opt_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the geometry optimization.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """

    return False


def _has_irc_nonconvergence_error_message(output_str):
    """ Assess whether the output file string contains the
        message signaling the failure of the
        Intrinsic Reaction Coordinate search.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: bool
    """


    return False 


ERROR_READER_DCT = {
    elstruct.par.Error.SCF_NOCONV: _has_scf_nonconvergence_error_message,
    elstruct.par.Error.MCSCF_NOCONV: lambda _: False,
    elstruct.par.Error.CC_NOCONV: lambda _: False,
    elstruct.par.Error.OPT_NOCONV: _has_opt_nonconvergence_error_message,
    elstruct.par.Error.IRC_NOCONV: _has_irc_nonconvergence_error_message,
    elstruct.par.Error.LIN_DEP_BASIS: lambda _: False
}

SUCCESS_READER_DCT = {
    elstruct.par.Success.SCF_CONV: _has_scf_convergence_message,
    elstruct.par.Success.OPT_CONV: _has_opt_convergence_message,
    elstruct.par.Success.IRC_CONV: _has_irc_convergence_message
}


def error_list():
    """ Constructs a list of errors that be identified from the output file.
    """
    return tuple(sorted(ERROR_READER_DCT.keys()))


def success_list():
    """ Constructs a list of successes that be identified from the output file.
    """
    return tuple(sorted(SUCCESS_READER_DCT.keys()))


def has_error_message(error, output_dct):
    """ Assess whether the output dictionary contains error messages
        for any of the procedures in the job.

        :param error: The type of error to check for
        :type error: elstruct.par.Error
        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    assert error in error_list()

    error_reader = ERROR_READER_DCT[error]
    if isinstance(error_reader, bool):
        err_val = False
    else:
        err_val = error_reader(output_dct)

    return err_val


def check_convergence_messages(error, success, output_dct):
    """ Assess whether the output dictionary indicates
        successful convergence of all requested procedures.

        :param error: The type of error to check for
        :type error: elstruct.par.Error
        :param success: The type of success to check for
        :type success: elstruct.par.Success
        :param output_dct: Dictionary of calculation results
        :type output_dct: dict
        :rtype: bool
    """
    assert error in error_list()
    assert success in success_list()

    if has_error_message(error, output_dct):
        job_success = SUCCESS_READER_DCT[success](output_dct)
    else:
        job_success = True

    return job_success