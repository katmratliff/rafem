RAFEM
=====

Morphodynamic river avulsion module designed to be coupled with
the Coastline Evolution Model (CEM) and Sedflux 3D.

The most current development version available from our git
repository:

http://github.com/katmratliff/avulsion-bmi


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
