====================================================
RAFEM: River Avulsion and Floodplain Evolution Model
====================================================

.. image:: https://img.shields.io/travis/sequence-dev/sequence.svg
        :target: https://travis-ci.org/sequence-dev/sequence

.. image:: https://ci.appveyor.com/api/projects/status/380ox1dv8hekefq9?svg=true
    :target: https://ci.appveyor.com/project/mcflugen/sequence/branch/develop

.. image:: https://readthedocs.org/projects/sequence/badge/?version=latest
        :target: https://sequence.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

About
-----

The River Avulsion and Floodplain Evolution Model (RAFEM) is a morphodynamic
model designed to be coupled with the Coastline Evolution Model (CEM).

Documentation is currently being updated (5/9/19). Please see our paper,
`Exploring Wave and Sea‐Level Rise Effects on Delta Morphodynamics With a Coupled River‐Ocean Model <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JF004757>`_ for a detailed model description. 


Requirements
------------

*Rafem* requires Python 3.

Apart from Python, *Rafem* has a number of other requirements, all of which
can be obtained through either *pip* or *conda*, that will be automatically
installed when you install *Rafem*.

To see a full listing of the requirements, have a look at the project's
*requirements.txt* file.

If you are a developer of *Rafem* you will also want to install
additional dependencies for running *Rafem*'s tests to make sure
that things are working as they should. These dependencies are listed
in *requirements-testing.txt*.

Installation
------------

To install *Rafem*, first create a new environment in
which *Rafem* will be installed. This, although not necessary, will
isolate the installation so that there won't be conflicts with your
base *Python* installation. This can be done with *conda* as::

  $ conda create -n rafem python=3
  $ conda activate rafem

Stable Release
++++++++++++++

*Rafem*, and its dependencies, can be installed either with *pip*
or *conda*. Using *pip*::

    $ pip install rafem

Using *conda*::

    $ conda install rafem -c conda-forge

From Source
+++++++++++

After downloading the *Rafem* source code, run the following from
*Rafem*'s top-level folder (the one that contains *setup.py*) to
install *Rafem* into the current environment::

  $ pip install -e .

Input Files
-----------

Rafem Parameter File
++++++++++++++++++++

The main *Rafem* input file is a yaml-formatted text file that lists
parameter values for the various components. Running the following will
print a sample *Rafem* parameter file::

  $ rafem show rafem

..code :: yaml

  shape:
  - 120
  - 200
  spacing:
  - 0.1
  - 0.1
  n0: 5.0
  nslope: 0.001
  max_rand: 0.1
  days: 7
  dt_day: 0.01
  rand_seed: 623
  Initial_SL: 0.0
  SLRR_m: 0.0
  SubRate_m: 0.0
  Sub_Start: 0
  ch_width: 10.0
  ch_depth: 1.0
  ch_discharge: 10.0
  A: 1.0
  c_f: 0.01
  C_0: 1.0
  sed_sg: 2.65
  init_cut_frac: 1
  super_ratio: 1.0
  short_path: 1
  WL_Z: 0.0
  WL_dist: 0
  blanket_rate_m: 0.0
  fine_dep_frac: 0.0
  splay_type: 2
  saveavulsions: false
  savecourseupdates: false

Output Files
------------

There are three main sets of output files. These are writen to the 
*output* folder as the model is running.
*  *output/elevation*: elevations of the entire model grid.
*  *output/profile*: elevations along the river profile
*  *output/river*: x, and y coordinates of the river profile

Each of these files is a CSV formatted text file. To create a plot
of one of these output files, use the *plot* subcommand. For example::

  $ rafem plot elevation

will plot the final elevations for the simulation in the current directory.
Use *rafem plot --help* to see further options.

Examples
--------

To run a simulation using the sample input files described above, you first
need to create a set of sample files. This can be done by hand or by running
`rafem setup` to get a default set of parameters that you can then edit.
For example::

  $ mkdir example
  $ cd example
  $ rafem setup

This command has created a new file, *rafem.yaml*, that you can edit for your
simulation.  To run *rafem* using this file::

  $ rafem run

This will have create a new folder, *output*, that contains the output files.
You can look at some of the output with the *plot* subcommand. For example,
the following will create a plot the final elevations::

  $ rafem plot elevation

Use the *--help* option to get help about other command line options.
