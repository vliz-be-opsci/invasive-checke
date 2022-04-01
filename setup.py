#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Rory Meyer",
    author_email='rory.meyer@vliz.be',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Use WRIMS + MarineRegions to determine if a given organism+location is invasive or not.",
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='invasive_checker',
    name='invasive_checker',
    packages=find_packages(include=['invasive_checker', 'invasive_checker.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/rory/invasive_checker',
    version='0.1.0',
    zip_safe=False,
)
