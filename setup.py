#!/usr/bin/env python

from setuptools import find_packages
from numpy.distutils.core import setup, Extension

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'xarray',
    'numpy',
    'matplotlib',
    'pandas',
    'metpy',
    'cartopy',
    'parse',
]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="eurec4a",
    author_email='',
    classifiers=[
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python",
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Atmospheric Science'],
    description="Software for working with data from the BAS twin otter",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='',
    name='twinotter',
    packages=find_packages(include=['twinotter']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='',
    version='0.1',
)