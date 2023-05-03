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
    #"/home/ygeorgievski/build/crossrate/DD/mess mess.inp >> stdout.log &> stderr.log"
    #"/home/ygeorgievski/build/crossrate/mess mess.inp >> stdout.log &> stderr.log"
    "mess mess.inp >> stdout.log &> stderr.log"
)
MESSRATEV2 = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=8\n"
    "ulimit -c 0\n"
    #"/home/ygeorgievski/build/crossrate/mess_test mess.inp >> stdout.log &> stderr.log"
    "mess-v2 mess.inp >> stdout.log &> stderr.log"
)
VARECOF_MULTI = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "MPI=`which mpirun`\n"
    # 'MPI_OPTIONS="-machinefile machines"'
    # 'MPI_OPTIONS="-host b460"'
    'MPI_OPTIONS="-n {}"\n'
    "VARECOFEXE=/home/mulvihcr/automech/VaReCoF/build/multi\n\n"
    "$MPI $MPI_OPTIONS $VARECOFEXE tst.inp >& varecof.out"
)
VARECOF_MCFLUX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/mulvihcr/automech/VaReCoF/build/mc_flux "
    "mc_flux.inp >& mc_flux.out"
)
VARECOF_CONV_STRUCT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/mulvihcr/automech/VaReCoF/build/convert_struct "
    "tst.inp >& varecof_struct_conv.out"
)
VARECOF_CONV_MULTI = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/mulvihcr/automech/VaReCoF/build/convert_multi "
    "convert.inp >& varecof_multi_conv.out"
)
VARECOF_TSTCHECK = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/tst_check >& tst_check.out"
)
INTDER = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    f"{EXTERN_PATH}/INTDER < intder.inp >& intder.out"
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
QCHEM5 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "export QC=/soft/qchem/5.3.0_mpich3/\n"
    "export QCSCRATCH=/scratch/$USER\n"
    "qchem -nt {} run.inp run.out"
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
    # INTDER
    'intder': INTDER,
    # VaReCoF
    'varecof_multi': VARECOF_MULTI,
    'varecof_conv_struct': VARECOF_CONV_STRUCT,
    'varecof_conv_multi':  VARECOF_CONV_MULTI,
    'varecof_mcflux': VARECOF_MCFLUX,
    'varecof_tstchk': VARECOF_TSTCHECK,
    # Electronic Structure
    'gaussian09': G09,
    'qchem5': QCHEM5,
    'molpro2021': MOLPRO_2021,
    'molpro2021_mppx': MOLPRO_2021_MPPX,
    'molpro2015': MOLPRO,
    'molpro2015_mppx': MOLPRO_MPPX,
    'molpro2015_mr': MOLPRO_MREF,
    'molpro2015_mr_mppx': MOLPRO_MREF_MPPX,
    'psi4': PSI4
}
