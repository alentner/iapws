"""A setuptools based setup module for the if97 Python package

"""
import re
from setuptools import setup, find_packages

# get long description from README.rst
with open('README.md', mode='r') as readme:
    long_description = readme.read()

# get package metadata by parsing __meta__ module
with open('src/iapws/__meta__.py', mode='r') as source:
    content = source.read().strip()
    metadata = {key: re.search(key + r'\s*=\s*[\'"]([^\'"]*)[\'"]', content).group(1)
                for key in ['__pkgname__', '__version__', '__authors__', '__contact__',
                            '__license__', '__website__', '__description__']}

# core dependancies
DEPENDANCIES = ['numpy', 'scipy']
OPTIONAL_PKG = {'test': ['matplotlib', ], }

# Define the setup configuration
setup(
    name                 = metadata['__pkgname__'],
    version              = metadata['__version__'],
    author               = metadata['__authors__'],
    author_email         = metadata['__contact__'],
    description          = metadata['__description__'],
    license              = metadata['__license__'],
    keywords             = 'water steam properties IAPWS IF97',
    url                  = metadata['__website__'],
    packages             = find_packages(where='src'),
    package_dir          = {'': 'src'},
    include_package_data = True,
    long_description     = long_description,
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics :: Material Properties',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        ],
    entry_points         = { },
    install_requires     = DEPENDANCIES,
    extras_require       = OPTIONAL_PKG, 
    install_requirments = [ ],
    )
