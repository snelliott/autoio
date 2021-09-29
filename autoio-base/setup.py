""" Install Interfaces to MESS, CHEMKIN, VaReCoF, ProjRot, and ThermP
"""

from distutils.core import setup

setup(
    name="autoio-base",
    version="0.8.0",
    packages=['autoparse',
              'ioformat',
              'autoread',
              'autoread._zmat',
              'autowrite'],
    package_dir={
        'autoparse': 'autoparse',
        'ioformat': 'ioformat',
        'autoread': 'autoread',
        'autowrite': 'autowrite'}
)
