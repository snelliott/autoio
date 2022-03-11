""" Library of BASH submission scripts for various programs
"""

import os

PATH = os.path.dirname(os.path.realpath(__file__))
EXTERN_PATH = os.path.join(PATH, 'extern')


PROJROT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "RPHt.exe >& rpht.out"
)
MESSPF = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=10\n"
    "ulimit -c 0\n"
    "messpf pf.inp >& pf.out"
)
MESSRATEV1 = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=8\n"
    "ulimit -c 0\n"
    "mess mess.inp >> stdout.log &> stderr.log"
)
MESSRATEV2 = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=8\n"
    "ulimit -c 0\n"
    "mess-v2 mess.inp >> stdout.log &> stderr.log"
)
VARECOF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "MPI=`which mpirun`\n"
    # 'MPI_OPTIONS="-machinefile machines"'
    # 'MPI_OPTIONS="-host b460"'
    # 'MPI_OPTIONS="-n {}"\n'
    'MPI_OPTIONS="-n 5"\n'
    "VARECOFEXE=/lcrc/project/CMRP/amech/VaReCoF/build/multi\n\n"
    "$MPI $MPI_OPTIONS $VARECOFEXE tst.inp >& varecof.out"
)
INTDER = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    f"{EXTERN_PATH}/INTDER < intder.inp >& intder.out"
)
MCFLUX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/lcrc/project/CMRP/amech/VaReCoF/build/mc_flux "
    "mc_flux.inp >& mc_flux.out"
)
VARECOF_CONV_STRUCT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/lcrc/project/CMRP/amech/VaReCoF/build/convert_struct "
    "tst.inp >& varecof_conv.out"
)
TSTCHECK = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/tst_check >& tst_check.out"
)
THERMP = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "thermp"
)
PAC99 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "pac99 << EOF >& pacc.out\n"
    "{}\n"
    "EOF"
)
G09 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "g09 run.inp run.out >> stdout.log &> stderr.log"
)
PSI4 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "psi4 -i run.inp -o run.out -n 8 >> stdout.log &> stderr.log"
)
MOLPRO_2021 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgievski/molpro_2021.2/bin/molpro "
    "-n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_2021_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgievski/molpro_2021.2/bin/molpro "
    "--mppx -n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n {} run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)

SCRIPT_DCT = {
    # MESS
    'messpf': MESSPF,
    'messrate-v1': MESSRATEV1,
    'messrate-v2': MESSRATEV2,
    # PAC99
    'pac99': PAC99,
    # ProjRot
    'projrot': PROJROT,
    # ThermP
    'thermp': THERMP,
    # VaReCoF
    'varecof': VARECOF,
    'varecof_conv_struct': VARECOF_CONV_STRUCT,
    'intder': INTDER,
    'mcflux': MCFLUX,
    'tstchk': TSTCHECK,
    # Electronic Structure
    'gaussian09': G09,
    'molpro2021': MOLPRO_2021,
    'molpro2021_mppx': MOLPRO_2021_MPPX,
    'molpro2015': MOLPRO,
    'molpro2015_mppx': MOLPRO_MPPX,
    'molpro2015_mr': MOLPRO_MREF,
    'molpro2015_mr_mppx': MOLPRO_MREF_MPPX,
    'psi4': PSI4
}
