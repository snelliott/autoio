".new_rates.py.swp""""
  Reads the output of a MESS calculation for the
  high-pressure and pressure-dependent rate constants
  corresponding to a given reaction.
"""

import sys
import numpy
import pandas as pd
import copy
from phydat import phycon
import autoparse.find as apf
from mess_io.reader._label import relabel
from mess_io.reader._label import name_label_dct
from autoreact.ktp_xarray import xarray_wrappers

# Global lists
UNWANTED_RXN_TYPS = ('fake', 'self', 'loss', 'capture', 'reverse')
DIRECTION_DCT = {'temp': 1000.0, 'pressure': 1.0, 'thresh': 1e-14}

# Functions for getting k(T,P) values from main MESS `RateOut` file
def get_rxn_ktp_dct(out_str,
                    filter_kts=True,
                    filter_reaction_types=UNWANTED_RXN_TYPS,
                    relabel_reactions=True,
                    tmin=None, tmax=None,
                    pmin=None, pmax=None, convert=True):
    """ Read a ktp dictionary for each reaction listed in the MESS output
        file string that requested by the user.

        Pressures in atm.
        K(T)s in cm3/mol.s [bimol] or 1/s [unimol]

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :param label_dct: dictionary to map name from MESS name to alternative
        :type label_dct: dict[str: str]
        :param filter_kts: filter unphysical, insignificant rate constants
        :type filter_kts: bool

        :rtype dict[float: (float, float)]
    """

    # Get the MESS rxn in the tuple format ((rct,), (prd,), third_body))
    # For each rxn pair, get rate constants, with filtering as indicated
    rxns = reactions(out_str)
    rxn_ktp_dct = {}
    for rxn in rxns:
        rxn_ktp_dct[rxn] = get_ktp(out_str, rxn[0][0], rxn[1][0],
                                   filter_kts=filter_kts, tmin=tmin,
                                   tmax=tmax, pmin=pmin, pmax=pmax,
                                   convert=convert)

    # Read the reactions; filter them if requested
    if filter_reaction_types:
        rxn_ktp_dct = filter_rxn_ktp_dct(
            rxn_ktp_dct,
            filter_fake=('fake' in filter_reaction_types),
            filter_self=('self' in filter_reaction_types),
            filter_loss=('loss' in filter_reaction_types),
            filter_capture=('capture' in filter_reaction_types),
            filter_reverse=('reverse' in filter_reaction_types)
        )

    # Reformat the dictionary keys to follow the tuple of tuples format
    if relabel_reactions:
        lbl_dct = name_label_dct(out_str)
        if lbl_dct is not None:
            rxn_ktp_dct = relabel(rxn_ktp_dct, lbl_dct)
        else:
            rxn_ktp_dct_relabeled = {}
            for key, val in rxn_ktp_dct.items():
                new_key = (tuple(key[0][0].split('+')),
                           tuple(key[1][0].split('+')), key[2])
                rxn_ktp_dct_relabeled[new_key] = val
            rxn_ktp_dct = copy.deepcopy(rxn_ktp_dct_relabeled)

    return rxn_ktp_dct


def get_ktp(output_str, reactant, product, filter_kts=True, tmin=None,
            tmax=None, pmin=None, pmax=None, convert=True):
    """ Parses the MESS output file string for the rate constants [k(T)]s
        for a single reaction for rate constants at all computed pressures,
        including the high-pressure limit.

        Pressures in atm.
        K(T)s in cm3/mol.s [bimol] or 1/s [unimol]

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :param reactant: label for the reactant used in the MESS output
        :type reactant: str
        :param product: label for the product used in the MESS output
        :type product: str
        :rtype dict[float: (float, float)]
    """

    # Get the MESS output lines
    out_lines = output_str.splitlines()
    # Initialize dictionary with high-pressure rate constants
    _highp = _highp_kts(out_lines, reactant, product)

    if _highp is not None:
        ktp = _highp
    else:
        ktp = {}

    # Read the pressures and convert them to atm if needed
    press, _ = get_pressures(output_str, mess_file='out')
    press = list(press)
    press.insert(0, numpy.inf)
    press.pop(len(press)-1)
    press = tuple(press)
    temps, _ = get_temperatures(output_str, mess_file='out')
    temporary_ktp = numpy.ndarray((len(press), len(temps)))

    temporary_ktp[0] = ktp

    # Update the dictionary with the pressure-dependent rate constants
    for pres_idx, pres in enumerate(_pres for _pres in press if not numpy.isinf(_pres)):
        temporary_ktp[pres_idx + 1] = _pdep_kts(out_lines, reactant, product, pres)
    bimol = (reactant[0] == 'P') or ('+' in reactant)

    full_ktp_dct = xarray_wrappers.from_data(list(temps), list(press), temporary_ktp)
    print(full_ktp_dct)

    # NVJY Note: convert_units not fully fixed to replace usage of ktp_dcts!
    # Currently commented so as to allow new_test_read_rates to function.

    #Note: filtering is before unit conversion, so bimolthresh is in cm^3.s^-1

    #full_ktp_dct = xarray_wrappers.dict_from_xarray(full_ktp_dct)

    if filter_kts:
        full_ktp_dct = filter_ktp(full_ktp_dct, bimol, tmin=tmin, tmax=tmax,
                               pmin=pmin, pmax=pmax)
    if convert:
        full_ktp_dct = convert_units(full_ktp_dct, bimol)

    return full_ktp_dct


def _highp_kts(out_lines, reactant, product):
    """ Parses the MESS output file string for the rate constants [k(T)]s
        for a single reaction at the high-pressure limit.

        :param out_lines: all of the lines of MESS output
        :type out_lines: list(str)
        :param reactant: label for the reactant used in the MESS output
        :type reactant: str
        :param product: label for the product used in the MESS output
        :type product: str
        :return rate_constants: all high-P rate constants for the reaction
        :rtype list(float)
    """

    # Find where the block of text where the high-pressure rates exist
    block_str = ('High Pressure Rate Coefficients ' +
                 '(Temperature-Species Rate Tables):')
    for i, line in enumerate(out_lines):
        if block_str in line:
            block_start = i
            for j in range(i, len(out_lines)):
                if '_________________________________' in out_lines[j]:
                    block_end = j
                    break
            break
    block_lines = out_lines[block_start: block_end]
    # Get the high-pressure rate constants
    rate_constants = None
    for i, line in enumerate(block_lines):
        # Find line with reactant name
        # If they match desired reactant, read rates
        # Rates start one line later in output
        if 'Reactant =' in line:
            mess_reac = line.strip().split()[2]
            if mess_reac == reactant:
                rate_constant = _parse_reactant_rate_constants(
                    block_lines, i+1, product)
                break
    
    return list(rate_constant)


def _pdep_kts(out_lines, reactant, product, pressure):
    """ Parses the MESS output file string for the rate constants [k(T,P)]s
        for a single reaction at a given numerical pressure, P.

        :param out_lines: all of the lines of MESS output
        :type out_lines: list(str)
        :param reactant: reactant name in MESS output
        :type reactant: str
        :param product: product name in MESS output
        :type product: str
        :param pressure: pressure that k(T,P)s will be read for
        :type pressure: float
        :return rate_constants: k(T,P)s for the reaction at given pressure
        :rtype list(float)
    """

    for i, line in enumerate(out_lines):
        if 'Temperature-Species Rate Tables:' in line:
            block_start = i
            for j in range(i, len(out_lines)):
                if '_________________________________' in out_lines[j]:
                    block_end = j
                    break
            break
    block_lines = out_lines[block_start: block_end]

    # Find where the block of text where the pressure-dependent rates exist
    pdep_vals = []
    for i, line in enumerate(block_lines):
        # Find line with reactant name and pressure
        # If they match desired reactant and pressure, read rates
        # Rates start two lines later in output
        if 'Reactant =' in line:
            tmp = line.strip().split()
            mess_reac = tmp[2]
            mess_press, mess_punit = float(tmp[5]), tmp[6]
            if (
                mess_reac == reactant and
                numpy.isclose(mess_press, pressure)
            ):
                atm_pressure = _convert_pressure(pressure, mess_punit)
                pdep_val = _parse_reactant_rate_constants(
                    block_lines, i+2, product)
                break

    return list(pdep_val)


def _parse_reactant_rate_constants(out_lines, block_start, product):
    """ Parses specific rate constants from the correct column
        in the MESS output file string.

        :param out_lines: all of the lines of MESS output
        :type out_lines: list(str)
        :param block_start: line num corresponding to reaction and pressure
        :type block_start: int
        :param product: label for the product used in the MESS output
        :type product: str
        :return rate_constants: all rate constants for the reaction
        :rtype: list(str, float)
    """

    # Find the column corresponding to the reaction
    product_col = 0
    product_headers = out_lines[block_start].strip().split()
    for i, product_header in enumerate(product_headers):
        if product == product_header:
            product_col = i
            break

    # Parse the following lines and store the constants in a list
    kts = []
    for i in range(block_start+1, len(out_lines)):
        if out_lines[i].strip() == '':
            break
        kts.append(out_lines[i].strip().split()[product_col])

    # Convert temps and rate constants to floats and combine values
    # only do so if the rate constant is defined (i.e., not '***')
    fin_kts = ()
    for kt_i in kts:
        new_kt = float(kt_i) if kt_i != '***' else numpy.nan
        fin_kts += (new_kt,)

    return fin_kts


# Functions for getting k(E)s and density-of-states from
# main MESS `MicroRateOut` file
def ke_dct(output_str, reactant, product):
    """ Parses the MESS output file string for the microcanonical
        rate constants [k(E)]s for a single reaction.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :param reactant: label for the reactant used in the MESS output
        :type reactant: str
        :param product: label for the product used in the MESS output
        :type product: str
        :return rate_constants: all high-P rate constants for the reaction
        :rtype dict[float: float]
    """

    # Form string for reaction header using reactant and product name
    reaction = _reaction_header(reactant, product)

    # Break up the file into lines
    out_lines = output_str.split('\n')

    # Determine col idx where reaction is
    head = out_lines[0].replace('E, kcal/mol', 'E').replace('D, mol/kcal', 'D')
    headers = head.split()
    col_idx = headers.index(reaction)

    _ke_dct = {0.0: 0.0}
    for i, line in enumerate(out_lines):
        if i not in (0, 1):
            tmp = line.strip().split()
            if tmp:
                _ke_dct[float(tmp[0])] = float(tmp[col_idx])

    return _ke_dct


def dos_rovib(ke_ped_out, sp_labels='auto'):
    """ Read the microcanonical pedoutput file and extracts rovibrational density
        of states of each fragment as a function of the energy

        units: kcal/mol, and for dos mol/kcal

        :param ke_ped_out: string of lines of microcanonical rates output file
        :type ke_ped_out: str
        :param sp_labels: type of pedspecies labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output
                'auto' sets 'inp' if it finds lbl dct
        :type sp_labels: str
        :return dos_df: dataframe(columns:prod1, prod2, rows:energy [kcal/mol])
                        with the density of states
        :rtype dos_df: dataframe(float)
    """
    # get label dictionary
    lbl_dct = name_label_dct(ke_ped_out)
    if sp_labels == 'auto':
        sp_labels = 'out'*(not lbl_dct) + 'inp'*(not not lbl_dct)

    ke_lines = ke_ped_out.splitlines()

    i_in = apf.where_in(
        'Bimolecular fragments density of states, mol/kcal', ke_lines)[0]+2
    mess_labels = ke_lines[i_in-1].strip().split()[2:]
    en_dos_all = numpy.array(
        [line.strip().split() for line in ke_lines[i_in:]], dtype=float).T
    energy = en_dos_all[:][0]
    dos_all = en_dos_all[:][1:].T

    # relabel if necessary
    _labels = []
    if sp_labels == 'inp' or sp_labels == 'out':
        for sp in mess_labels:
            bim, frag_n = sp.split('_')
            if sp_labels == 'inp':
                sp_tosplit = lbl_dct[bim]
            elif sp_labels == 'out':
                sp_tosplit = bim
            try:
                _labels.append(sp_tosplit.split('+')[int(frag_n)])
            except IndexError:
                print('*Error: bimol species should be named as P1+P2 \
                    with P1, P2 being the fragment names')
                sys.exit()
    else:
        print('*Error: sp_labels must be "inp" (as in mess input) \
            or "out" (as in mess output)')
        sys.exit()

    dos_rovib_df = pd.DataFrame(dos_all, index=energy, columns=_labels)
    # drop potentially duplicate columns WARNING CHECK THE EFFECT OF THIS
    dos_rovib_df = dos_rovib_df.T.drop_duplicates().T

    return dos_rovib_df


# Functions for getting barrier heights and corresponding reactants, products
# not tested because currently unused
def energies(output_str):
    """ Parses the MESS output and gets the energies of all species
        :param output_str: string of lines of MESS output file
        :type output_str: str

        :return species_en, barriers_en: species and barrier energy in kcal/mol
        :rtype series(index=name), series(index=(reac,prod))
    """

    # Break up the file into lines
    out_lines = numpy.array(output_str.split('\n\n'))
    headers = ['Wells', 'Bimolecular Products',
               'Well-to-Bimolecular', 'Well-to-Well']
    # Determine where block is
    block_indexes = apf.where_in_any(headers, out_lines)

    # preallocations
    _species = []
    _species_ene = []
    _barriers = []
    _barriers_ene = []

    for i in block_indexes:
        block = out_lines[i]
        lines = block.splitlines()

        if 'Wells' in block or 'Bimolecular Products' in block:
            for line in lines[2:]:
                spc, ene = line.strip().split()[0:2]
                _species.append(spc)
                _species_ene.append(float(ene))

    species_ene_s = pd.Series(_species_ene, index=_species)

    for i in block_indexes:
        block = out_lines[i]
        lines = block.splitlines()
        if 'Well-to-Bimolecular' in block:
            for line in lines[2:]:
                _, ene, reac, _, prod = line.strip().split()
                _barriers.append((reac, prod))
                _barriers.append((prod, reac))
                # relative energy
                _barriers_ene.append(float(ene)-species_ene_s[reac])
                _barriers_ene.append(float(ene)-species_ene_s[prod])

        if 'Well-to-Well' in block:
            for line in lines[2:]:
                _, ene, reac, _, prod, _ = line.strip().split()
                _barriers.append((reac, prod))
                _barriers.append((prod, reac))
                _barriers_ene.append(float(ene)-species_ene_s[reac])
                _barriers_ene.append(float(ene)-species_ene_s[prod])

    barriers_ene_s = pd.Series(_barriers_ene, index=_barriers)

    return species_ene_s, barriers_ene_s


def barriers(barriers_ene_s, species_ene_s, reac, prod):
    """ Extract fw and bw energy barrier for reaction reac->prod
        if barrier not found, save the DE of the reaction
    """
    findreac = [reac == key[0] for key in barriers_ene_s.keys()]
    findprod = [prod == key[1] for key in barriers_ene_s.keys()]
    findreac_rev = [reac == key[1] for key in barriers_ene_s.keys()]
    findprod_rev = [prod == key[0] for key in barriers_ene_s.keys()]

    if (reac, prod) in barriers_ene_s.keys():
        dene_fw = barriers_ene_s[(reac, prod)]
        # DW_BW = barriers_ene_s[(prod, reac)]

    elif (
        any(findreac) or any(findreac_rev) and
        any(findprod) or any(findprod_rev)
    ):
        # check if you have reac->wr->wp->prod like in habs
        connect_reac = numpy.array(list(barriers_ene_s.keys()))[findreac]
        fromreac = [p[1] for p in connect_reac]
        connect_prod = numpy.array(list(barriers_ene_s.keys()))[findprod]
        fromprod = [r[0] for r in connect_prod]
        possible_index = [(p, r) for p in fromreac for r in fromprod]
        possible_index += [(reac, p) for p in fromprod]
        possible_index += [(r, prod) for r in fromreac]

        flag = 0

        for idx in possible_index:

            try:
                delta_ene = barriers_ene_s[idx]
                # barriers like r=>wr=>wp=>p
                dene_fw = delta_ene - barriers_ene_s[(idx[0], reac)]
                dene_bw = barriers_ene_s[(idx[1], idx[0])] - \
                    barriers_ene_s[(idx[1], prod)]
                flag = 1
            except (KeyError, ValueError):
                try:
                    # barriers like r=>wr=>p
                    # reacs=>wr
                    dene_fw = barriers_ene_s[(idx[0], prod)] - \
                        barriers_ene_s[(idx[0], reac)]

                    dene_bw = dene_fw + barriers_ene_s[(prod, idx[0])]
                    flag = 1
                except (KeyError, ValueError):
                    continue

        if flag == 0:
            print('*Warning: indirect connection between '
                  f'{reac} and {prod} not found')
            print('saving the DE instead \n')
            try:
                dene_fw = species_ene_s[prod] - species_ene_s[reac]
                dene_bw = species_ene_s[reac] - species_ene_s[prod]
            except KeyError:
                print('*Error: species '
                      f'{reac} and {prod} not found. Now exiting \n')
    else:
        print(f'*Error: species {reac} and {prod} not found')
        sys.exit()

    return dene_fw, dene_bw


# Functions to read temperatures and pressures
def get_temperatures(file_str, mess_file='out'):
    """ Get temps
    """

    if mess_file == 'out':
        temps = _temperatures_output_string(file_str)
    elif mess_file == 'inp':
        temps = _temperatures_input_string(file_str)
    else:
        temps = ()

    return list(temps)


def get_pressures(file_str, mess_file='out'):
    """ Get pressures
    """

    if mess_file == 'out':
        _pressures = _pressures_output_string(file_str)
    elif mess_file == 'inp':
        _pressures = _pressures_input_string(file_str)
    else:
        _pressures = ()

    return list(_pressures)


def _temperatures_input_string(input_str):
    """ Reads the temperatures from the MESS input file string
        that were used in the master-equation calculation.
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return temperatures: temperatures in the input
        :rtype: list(float)
        :return temperature_unit: unit of the temperatures in the input
        :rtype: str
    """

    # Get the MESS input lines
    mess_lines = input_str.splitlines()
    for line in mess_lines:
        if 'TemperatureList' in line:
            temps = tuple(float(val) for val in line.strip().split()[1:])
            temp_unit = line.strip().split('[')[1].split(']')[0]
            break

    return temps, temp_unit


def _temperatures_output_string(output_str):
    """ Reads the temperatures from the MESS output file string
        that were used in the master-equation calculation.
        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return temperatures: temperatures in the output
        :rtype: list(float)
        :return temperature_unit: unit of the temperatures in the output
        :rtype: str
    """

    # Get the MESS output lines
    mess_lines = output_str.splitlines()

    # Find the block of lines where the temperatures can be read
    for i, line in enumerate(mess_lines):
        if 'Pressure-Species Rate Tables:' in line:
            block_start = i
            for j in range(i, len(mess_lines)):
                if 'Temperature-Species Rate Tables:' in mess_lines[j]:
                    block_end = j
                    break
            break

    # Read the temperatures
    temps = []
    for i in range(block_start, block_end):
        if 'Temperature =' in mess_lines[i]:
            tmp = mess_lines[i].strip().split()
            temps.append(float(tmp[5]))
    temps = list(set(temps))
    temps.sort()

    # Read unit
    for i in range(block_start, block_end):
        if 'Temperature =' in mess_lines[i]:
            temp_unit = mess_lines[i].strip().split()[6]
            break

    return tuple(temps), temp_unit


def _pressures_input_string(input_str):
    """ Reads the pressures from the MESS input file string
        that were used in the master-equation calculation.
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return pressures: pressures in the input
        :rtype: list(str, float)
        :return pressure_unit: unit of the pressures in the input
        :rtype: str
    """

    # Get the MESS input lines
    mess_lines = input_str.splitlines()
    for line in mess_lines:
        if 'PressureList' in line:
            _pressures = [float(val) for val in line.strip().split()[1:]]
            pressure_unit = line.strip().split('[')[1].split(']')[0]
            break

    # Append high pressure
    _pressures.append('high')

    return tuple(_pressures), pressure_unit


def _pressures_output_string(output_str):
    """ Reads the pressures from the MESS output file string
        that were used in the master-equation calculation.
        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return pressures: pressures in the output
        :rtype: list(str, float)
        :return pressure_unit: unit of the pressures in the output
        :rtype: str
    """

    # Get the MESS output lines
    mess_lines = output_str.splitlines()

    # Find the block of lines where the pressures can be read
    pressure_str = 'Pressure-Species Rate Tables:'
    for i, line in enumerate(mess_lines):
        if pressure_str in line:
            block_start = i
            break

    # Read the pressures
    _pressures = []
    for i in range(block_start, len(mess_lines)):
        if 'P(' in mess_lines[i]:
            pressure_unit = mess_lines[i].strip().split('(')[1].split(')')[0]
            pressure_start = i+1
            for j in range(pressure_start, len(mess_lines)):
                if 'O-O' in mess_lines[j]:
                    break
                tmp = mess_lines[j].strip().split()
                _pressures.append(float(tmp[0]))
            break

    # Append high pressure
    _pressures.append(numpy.inf)

    return tuple(_pressures), pressure_unit


def _convert_pressure(pressure, pressure_unit):
    """ Convert a set of pressures using the pressure unit, accounting
        for high-pressure in list

        :param pressure: pressures to convert
        :type pressure: tuple(float, str)
        :param pressure_unit: unit of pressures
        :type pressure_unit: str
        :rtype: tuple(float, str)
    """
    
    if pressure != 'high':
        if pressure_unit == 'atm':
            conv = 1.0
        elif pressure_unit == 'torr':
            conv = phycon.TORR2ATM
        elif pressure_unit == 'bar':
            conv = phycon.BAR2ATM
        pressure *= conv

    return pressure


# Read the reactions
def reactions(out_str, third_body=(None,)):
    """ Read the reactions from the output file.

        Currently, we assume we don't care about the third body
        and have the user pass it.

        Returns the reactions in tuple list (reacs, prods, third_body)
        where reacs and prods are tuples

        Ignores 'Capture' reactions
    """

    # Split into lines
    out_lines = out_str.splitlines()

    # Grab reactions from T-dep
    for i, line in enumerate(out_lines):
        if 'Temperature-Species Rate Tables:' in line:
            block_start = i
            for j in range(i, len(out_lines)):
                if 'Temperature-Pressure Rate Tables' in out_lines[j]:
                    block_end = j
            break
    reac_lines = out_lines[block_start:block_end+1]

    # Read all of the reactions out of the file
    # Reactant = <reactant>
    #
    # T(K) <prod1> <prod2> ... Loss Capture
    rxns = ()
    for i, line in enumerate(reac_lines):
        if 'Reactant = ' in line and 'Pressure =' in line:
            _reac = line.strip().split()[2]
            _prods = reac_lines[i+2].strip().split()[1:]
            for _prod in _prods:
                rxns += (
                    ((_reac,), (_prod,), third_body),
                )

    # Remove duplcates while preserving order
    rxns = tuple(n for i, n in enumerate(rxns)
                 if n not in rxns[:i])

    return rxns


def filter_rxn_ktp_dct(rxn_ktp_dct,
                       filter_fake=True,
                       filter_self=True,
                       filter_loss=True,
                       filter_capture=True,
                       filter_reverse=True):
    """ Filter the reactions from a ktp dictionary
    """

    filt_rxn_ktp_dct = {}
    for rxn, ktp_dct in rxn_ktp_dct.items():
        rct, prd, tbody = rxn
        if prd == ('Loss',):
            if filter_loss:
                continue
        if prd == ('Capture',):
            if filter_capture:
                continue
        # if (  # can we delete this? 
        #     any('F' in rgt for rgt in rct+prd)  or
        #     any('FW' in rgt for rgt in rct+prd)
        # ):
        if any('Fake' in rgt for rgt in rct+prd):
            if filter_fake:
                continue
        if rct == prd:
            if filter_self:
                continue
        if filter_reverse:
            is_desired = is_desired_direction(rxn, rxn_ktp_dct)
            if not is_desired:
                continue
        # If continues not hit, reaction good to be added to new dct
        filt_rxn_ktp_dct[rxn] = ktp_dct

    return filt_rxn_ktp_dct


# Helper functions
def _reaction_header(reactant, product):
    return reactant + '->' + product


def is_desired_direction(rxn, rxn_ktp_dct, direction_dct=DIRECTION_DCT):
    """ Decides whether the current direction is the desired one by comparing
        forward and backward rate constants

        :rtype: bool
    """

    def get_avg(ktp_dct, targ_temp, pressures):
        """ Gets a singular rate constant at a desired temp and pressure
        """
        rate_sum = 0
        nrates = 0  # number of forward rates (i.e., not NaNs)
        # Get forward rates
        for pressure in pressures:
            curr_temps = ktp_dct[pressure][0]
            tidx = (numpy.abs(numpy.array(curr_temps) - targ_temp)).argmin()
            curr_rate = ktp_dct[pressure][1][tidx]
            if not numpy.isnan(curr_rate):
                rate_sum += curr_rate
                nrates += 1
        # Average
        if nrates == 0:  # if nothing was found!
            rate_avg = 1e-100
        else:
            rate_avg = rate_sum / nrates

        return rate_avg

    # Load things
    rct, prd, tbody = rxn
    targ_temp = direction_dct['temp']
    thresh = direction_dct['thresh']
    fwd_ktp_dct = rxn_ktp_dct[rxn]
    bck_ktp_dct = rxn_ktp_dct[(prd, rct, tbody)]
    if fwd_ktp_dct == {} and bck_ktp_dct == {}:
        print(f'Error: both directions empty for the rxn {rxn}')
    #if fwd_ktp_dct == {}:
    #    return False
    #if bck_ktp_dct == {}:
    #    return True

    # Determine the molecularity of the reaction
    unimol_unimol = False
    bimol_bimol = False
    unimol_bimol = False
    if 'W' in rct[0] and 'W' in prd[0]:  # unimol to unimol
        unimol_unimol = True
    if 'P' in rct[0] and 'P' in prd[0]:  # bimol to bimol
        bimol_bimol = True

    # Get rate value to use for comparison
    # If both directions have 'high' values, use that at target T
    if 'high' in fwd_ktp_dct.keys() and 'high' in bck_ktp_dct.keys():
        fwd_avg = get_avg(fwd_ktp_dct, targ_temp, ['high'])
        bck_avg = get_avg(bck_ktp_dct, targ_temp, ['high'])
    # Otherwise, average over all pressures except 'high' at target T
    else: 
        fwd_pressures = [key for key in fwd_ktp_dct.keys() if key != 'high']
        bck_pressures = [key for key in bck_ktp_dct.keys() if key != 'high']
        fwd_avg = get_avg(fwd_ktp_dct, targ_temp, fwd_pressures)
        bck_avg = get_avg(bck_ktp_dct, targ_temp, bck_pressures)

    # Finally, determine whether the current direction is desired
    # If unimol > unimol or bimol > bimol, larger direction is preferred 
    is_desired = False
    if unimol_unimol or bimol_bimol:
        if fwd_avg > bck_avg:
            is_desired = True
    # If bimol > unimol (or reverse), consider the threshold value
    else:
        if 'W' in rct[0]:  # written as unimol to bimol
            if bck_avg / phycon.NAVO < thresh:  # if bimol less than thresh
                    is_desired = True
        else:  # written as bimol to unimol
            if fwd_avg / phycon.NAVO > thresh:  # if bimol greater than thresh
                is_desired = True
        
    return is_desired
            

def filter_ktp(ktp, bimol,
                   tmin=None, tmax=None, pmin=None, pmax=None):
    """ Filters out bad or undesired rate constants from a ktp dictionary
    """

    def get_valid_k(kts, bimol, tmin = None, tmax = None, 
                    bimolthresh=1.0e-24):
        """ Takes in lists of temperature-rate constant pairs [T,k(T)]
            and removes invalid pairs for which
            (1) k(T) < 0
            (3) k(T) < some threshold for bimolecular reactions, or
            (4) T is outside the cutoff
            (5) there are more than 2 oscillations from pos. to neg.

            :param temps: temperatures at which rate constants are defined (K)
            :type temps: list(float)
            :param kts: rate constants (s-1 or cm^3.s-1)
            :type kts: list(str, float)
            :param bimol: whether or not the reaction is bimolecular
            :type bimol: Bool
            :param bimolthresh: threshold below which bimolecular rate
                constants will be removed (cm^3.s^-1)
            :type bimolthresh: float
            :return valid_t: List of vaild temperatures
            :return valid_k: List of vaild rate constants
            :rtype: (numpy.ndarray, numpy.ndarray)
            """

        bimolthresh = 1.0e-24

        # Set min temperature to user input, if none use either
        # min of input temperatures or
        # if negative kts are found, set min temp to be just above highest neg.
        #imax_neg_idx = None
        #sing_kts = []  # empty list that won't have any Nones
        #for kt_idx, sing_kt in enumerate(kts):
        #    # find idx for max temperature for which kt is negative, if any
        #    float_sing_kt = float(sing_kt) if sing_kt is not None else 1.0
        #    if float_sing_kt < 0.0:
        #        max_neg_idx = kt_idx
        #    sing_kts.append(float_sing_kt) 

        #   Otherwise, use requested tmin value or minimum of input temps
        # Set tmin to None (i.e., no valid k(T)) if highest T has neg. k(T)
        
        if tmin is not None and tmax is not None:
            # Check for sign changes
            allowed = 2  # number of allowed changes
            nsign_chgs = numpy.count_nonzero((numpy.diff(numpy.sign(kts)) != 0) * 1)
            if nsign_chgs > allowed:
                valid_k = []

        # Grab the temperature, rate constant pairs which correspond to
        # temp > 0, temp within tmin and tmax, rate constant defined (not ***)
        kthresh = 0.0 if not bimol else bimolthresh
        kts[kts < kthresh] = numpy.nan

        return kts

    # Filter the kts based on temps, negatives, None, and bimolthresh
    if tmin == None:
        tmin = float(ktp.coords['temp'].min(skipna=True))
    ktp = ktp.where(ktp['temp'] >= tmin, numpy.nan)
    if tmax == None:
        tmax = float(ktp.coords['temp'].max(skipna=True))
    ktp = ktp.where(ktp['temp'] <= tmax, numpy.nan)

    temp_idx = 0
    for temp in xarray_wrappers.get_temperatures(ktp):
        kts = xarray_wrappers.get_tslice(ktp, temp)
        kts = kts.to_numpy()
        filt_kts = get_valid_k(kts, bimol)
        ktp[:, temp_idx] = filt_kts
        temp_idx += 1

    # Remove undesired pressures if pmin and/or pmax were given (leaves 'high'
    # untouched if it is present)
    if pmin == None:
        pmin = float(ktp.coords['pres'].min(skipna=True))
    ktp = ktp.where(ktp['pres'] >= pmin, numpy.nan)
    if pmax == None:
        pmax = float(ktp.pres[numpy.isfinite(ktp.pres)].max(skipna=True))
    ktp = ktp.where((ktp['pres'] <= pmax) | (ktp['pres'] == numpy.inf), numpy.nan)

    return ktp


def convert_units(ktp, bimol):
    """ Convert units from cm^3.s^-1 to cm^3.mol^-1.s^-1 if rxn is bimolecular
    """

    if bimol:
        ktp *= phycon.NAVO

    return ktp


def filter_reactions(rxns,
                     filter_fake=True,
                     filter_self=True,
                     filter_loss=True,
                     filter_capture=True,
                     filter_reverse=True):
    """ Filter the reactions from a ktp dictionary
        Leaving this function here for the sake of _wellextend.py
    """

    filt_rxns = ()
    for rxn in rxns:
        rct, prd, tbody = rxn
        if prd == ('Loss',):
            if filter_loss:
                continue
        if prd == ('Capture',):
            if filter_capture:
                continue
        # if (
        #     any('F' in rgt for rgt in rct+prd)  or
        #     any('FW' in rgt for rgt in rct+prd)
        # ):
        if any('Fake' in rgt for rgt in rct+prd):
            if filter_fake:
                continue
        if rct == prd:
            if filter_self:
                continue
        if filter_reverse:
            if (prd, rct, tbody) in filt_rxns:
                continue
        # If continues not hit, reaction good to be added to new dct
        filt_rxns += (rxn,)

    return filt_rxns
