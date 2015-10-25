#! /usr/bin/env python
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages


setup(name='bmi-avulsion',
      version='0.1.0',
      author='Katherine Ratliff',
      author_email='k.ratliff@duke.edu',
      description='BMI Avulsion Module',
      long_description=open('README.rst').read(),
      url='https://github.com/katmratliff/avulsion-bmi',
      license='MIT',
      packages=find_packages(),
)
