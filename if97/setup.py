"""A setuptools based setup module for the if97 Python package

"""
from setuptools import setup

# Get the long description from the README file
with open('README.rst') as file:

    long_description = file.read()

# Define the setup configuration
setup(
    name = 'if97',
    version = '1.0.4',
    description = 'A Python package for providing water properties', 
    long_description = long_description,
    url = 'https://github.com/alentner/if97',
    author = 'Aaron Lentner',
    author_email = 'aaron.d.lentner@gmail.com',
    license = 'MIT',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics :: Material Properties',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        ],
    keywords = 'water steam properties industrial formulation 97 IAPWS',
    packages=['if97'],
    install_requirments = [ ],
    entry_points = { },
    include_package_data = True,
    zip_safe = False
    )