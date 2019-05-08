#! /usr/bin/env python
from setuptools import find_packages, setup

import versioneer


setup(
    name="rafem",
    version=versioneer.get_version(),
    author="Katherine Ratliff",
    author_email="k.ratliff@duke.edu",
    description="River Avulsion Flooplain Evolution Model",
    long_description=open("README.rst").read(),
    url="https://github.com/katmratliff/avulsion-bmi",
    license="MIT",
    install_requires=[
        "bmipy",
        "click",
        "numpy",
        "matplotlib",
        "pyyaml",
        "scipy",
        "six",
    ],
    packages=find_packages(),
    entry_points={"console_scripts": ["rafem=rafem.main:main"]},
    cmdclass=versioneer.get_cmdclass(),
)
