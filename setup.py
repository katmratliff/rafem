#! /usr/bin/env python
from setuptools import find_packages, setup

import versioneer


setup(
    name="rafem",
    version=versioneer.get_version(),
    author="Katherine Ratliff",
    author_email="k.ratliff@duke.edu",
    description="River Avulsion Flooplain Evolution Model",
    long_description=open("README.rst", encoding="utf-8").read(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    url="https://github.com/katmratliff/avulsion-bmi",
    license="MIT",
    install_requires=open("requirements.txt", "r").read().splitlines(),
    packages=find_packages(),
    entry_points={"console_scripts": ["rafem=rafem.main:rafem"]},
    cmdclass=versioneer.get_cmdclass(),
)
