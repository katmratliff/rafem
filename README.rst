RAFEM
=====

The River Avulsion and Floodplain Evolution Model (RAFEM) is a morphodynamic model designed to be coupled with
the Coastline Evolution Model (CEM).

Documentation is currently being updated (5/9/19). Please see our paper, `Exploring Wave and Sea‐Level Rise Effects on Delta Morphodynamics With a Coupled River‐Ocean Model <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JF004757>`_ for a detailed model description. 


Install
-------

The sources for RAFEM can be downloaded from the Github repo.

You can either clone the public repository:

  .. code:: bash

    $ git clone git://github.com/katmratliff/rafem

Or download the tarball:

  .. code:: bash

    $ curl  -OL https://github.com/katmratliff/rafem/tarball/master

Once you have a copy of the source, you can install it with:

  .. code:: bash

    $ pip install -e .

Run
===

After installing RAFEM, you can run the model with the *rafem* command.

  .. code:: bash

    $ rafem input.yaml

Use the *--help* option to get help about other command line options.
