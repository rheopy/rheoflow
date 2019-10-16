# -*- coding: utf-8 -*-

# Learn more: https://github.com/rheopy/rheoflow

from setuptools import setup, find_packages

readme=''
licence='Apache 2'

setup(
    name='rheoflow',
    version='0.1.0',
    description='Python library for non-Newtonian flow calculations',
    long_description=readme,
    author='William Hartt',
    author_email='whartt@gmail.com',
    url='https://github.com/rheopy/rheoflow',
    license=license,
    install_requires=['numpy','matplotlib'],
    packages=['rheoflow']
)
