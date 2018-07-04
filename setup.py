#! /usr/bin/env python
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

from model_metadata.utils import get_cmdclass, get_entry_points


pymt_components = [
    (
        "BmiRiverModule=rafem:BmiRiverModule",
        ".bmi",
    )
]


setup(name='rafem',
      version='0.1.0',
      author='Katherine Ratliff',
      author_email='k.ratliff@duke.edu',
      description='River Avulsion Flooplain Evolution Model',
      long_description=open('README.rst').read(),
      url='https://github.com/katmratliff/avulsion-bmi',
      license='MIT',
      packages=find_packages(),
      cmdclass=get_cmdclass(pymt_components),
      entry_points=get_entry_points(pymt_components),
)
