"""
Writes the full sections of a MESS input file
"""

import os
from ioformat import build_mako_str
from ioformat import indent
from ioformat import remove_trail_whitespace
from phydat import phycon
from mess_io.writer._mol_inf import core_rigidrotor
from mess_io.writer._spc import molecule
from mess_io.writer._rxnchan import species


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')


# Write the full input file strings
def messrates_inp_str(globkey_str, rxn_chan_str,
                      energy_trans_str=None, well_lump_str=None):
    """ Combine various MESS strings together to combined MESS rates
    """

    # Create dictionary to fill template
    full_rxn_inp_keys = {
        'globkey_str': globkey_str,
        'energy_trans_str': energy_trans_str,
        'well_lump_str': well_lump_str,
        'rxn_chan_str': rxn_chan_str
    }

    mess_inp_str = build_mako_str(
        template_file_name='full_rates_inp.mako',
        template_src_path=SECTION_PATH,
        template_keys=full_rxn_inp_keys)

    return remove_trail_whitespace(mess_inp_str)


def messpf_inp_str(globkey_str, spc_str):
    """ Combine various MESS strings together to combined MESSPF
    """
    return '\n'.join([globkey_str, spc_str]) + '\n'


def messhr_inp_str(geo, hind_rot_str):
    """ Special MESS input string to calculate frequencies and ZPVEs
        for hindered rotors
    """

    global_pf_str = global_pf_input(
        temperatures=(100.0, 200.0, 300.0, 400.0, 500),
        rel_temp_inc=0.001,
        atom_dist_min=0.6)
    dat_str = molecule(
        core=core_rigidrotor(geo, 1.0),
        freqs=(1000.0,),
        elec_levels=((0.0, 1.0),),
        hind_rot=hind_rot_str,
    )
    spc_str = species(
        spc_label='Tmp',
        spc_data=dat_str,
        zero_ene=0.0
    )

    return messpf_inp_str(global_pf_str, spc_str)


# Write individual sections of the input file
def global_rates_input(temperatures, pressures,
                       excess_ene_temp=None, well_extend='auto'):
    """ Writes the global keywords section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param temperatures: List of temperatures (in K)
        :type temperatures: float
        :param pressures: List of pressures (in atm)
        :type pressures: float
        :return global_str: String for section
        :rtype: string
    """

    # Format temperature and pressure lists
    temperature_list = '  '.join(str(val) for val in temperatures)
    pressure_list = '  '.join(str(val) for val in pressures)

    # Format the other keywords as needed
    if excess_ene_temp is not None:
        assert isinstance(excess_ene_temp, float), (
            'ExcessEnergyOverTemperature value must be a float'
        )
        excess_ene_temp_str = f'{excess_ene_temp:.2f}'
    else:
        excess_ene_temp_str = None
    if well_extend is not None:
        if well_extend != 'auto':
            assert isinstance(well_extend, float), (
                'WellExtension value must be a float'
            )
            well_extend_str = f'{well_extend:.2f}'
        else:
            well_extend_str = ''
    else:
        well_extend_str = None

    # Create dictionary to fill template
    globrxn_keys = {
        'temperatures': temperature_list,
        'pressures': pressure_list,
        'excess_ene_temp': excess_ene_temp_str,
        'well_extend': well_extend_str,
    }

    return build_mako_str(
        template_file_name='global_rates.mako',
        template_src_path=SECTION_PATH,
        template_keys=globrxn_keys)


def global_pf_input(temperatures=(),
                    temp_step=100, ntemps=30,
                    rel_temp_inc=0.001, atom_dist_min=1.13384):
    """ Writes the global keywords section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param temperatures: List of temperatures (in K)
        :type temperatures: list(float)
        :param temp_step: temperature step (in K)
        :type temp_step: float
        :param ntemps: number of temperature values on grid
        :type ntemps: int
        :param rel_temp_inc: increment for temps
        :type rel_temp_inc: float
        :param atom_dist_min: cutoff for atom distances (Bohr)
        :type atom_dist_min: float
        :return global_pf_str: string for section
        :rtype: string
    """

    if temperatures:
        temperature_list = '  '.join(str(val) for val in temperatures)
        temp_step = None
        ntemps = None
    else:
        temperature_list = ''

    # Convert the atom distance minimum
    atom_dist_min = f'{atom_dist_min * phycon.BOHR2ANG:.2f}'

    # Create dictionary to fill template
    globpf_keys = {
        'temperatures': temperature_list,
        'temp_step': temp_step,
        'ntemps': ntemps,
        'rel_temp_inc': rel_temp_inc,
        'atom_dist_min': atom_dist_min
    }

    return build_mako_str(
        template_file_name='global_pf.mako',
        template_src_path=SECTION_PATH,
        template_keys=globpf_keys)


def global_energy_transfer_input(edown_str, collid_freq_str):
    """ Writes the global energy transfer section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param edown_str: String for the energy down parameters
        :type edown_str: str
        :param collid_freq_str: String for the collisional freq parameters
        :type collid_freq_str: str
        :rtype: str
    """

    edown_str = indent(edown_str, 2)
    collid_freq_str = indent(collid_freq_str, 2)

    # Create dictionary to fill template
    glob_etrans_keys = {
        'edown_str': edown_str,
        'collid_freq_str': collid_freq_str
    }

    return build_mako_str(
        template_file_name='global_etrans.mako',
        template_src_path=SECTION_PATH,
        template_keys=glob_etrans_keys)


# Write data to output file formats if you want to make a formatted output file
def pf_output(fml_str, temps, logq, dq_dt, d2q_dt2, svals=None, cpvals=None):
    """ Writes partition function data into a string that is formatted like the
        output file
    """

    mess_out_str = 'Natural log of the partition function '
    mess_out_str += 'and its derivatives:\n'
    mess_out_str += f' T, K            {fml_str}\n'
    for idx, _ in enumerate(temps):
        mess_out_str += '\n'
        mess_out_str += f'{temps[idx]:>8.3f}    '
        mess_out_str += f'{logq[idx]:>8.6f}    '
        mess_out_str += f'{dq_dt[idx]:>8.8f}    '
        mess_out_str += f'{d2q_dt2[idx]:>8.6e}    '
        if svals is not None:
            mess_out_str += f'{svals[idx]:>8.6f}    '
        if cpvals is not None:
            mess_out_str += f'{cpvals[idx]:>8.6f}    '

    return mess_out_str
