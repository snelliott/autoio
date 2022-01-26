""" Install AutoIO-Base library that wrap simple
    python I/O functions to perform more complex tasks
"""

from distutils.core import setup


setup(
    name="autoio-base",
    version="0.9.2",
    packages=[
        'autoparse',
        'ioformat',
        'autoread',
        'autoread._zmat',
        'autowrite'
    ],
    package_dir={
        'autoparse': 'autoparse',
        'ioformat': 'ioformat',
        'autoread': 'autoread',
        'autowrite': 'autowrite'
    }
)
