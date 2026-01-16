.. _install:

=======
Install
=======

To install sith, run these lines on your terminal:

.. code-block:: bash

    git clone git@github.com:graeter-group/sith.git
    cd sith
    pip install -e .

Note that "-e" flag allows to use sith in such a way that you can make
changes on it and use them immediately without having to reinstall the package
again. However, you can avoid this flag, which would generate a copy of sith in
the python packages path.

.. _load-modules:

------------
Load-modules
------------

So far, sith is based on gaussian because it has implemented the calculation of
forces directly in the DOFs. For future releases, we plan to expand the use of
sith to other QM softwares. We developed our code specifically for g09, but it
should work well with any newer gaussian version. Given that the gaussian
version is completely customizable, the implemented tools of sith that use
gaussian requiere a function that has to be defined in either of the two ways
shown bellow.

Also, in HPC (in clusters), it is common to load the necessary softwares as
modules. This step is also shown in the options bellow. If you do not need to
load modules, ignore that part.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Option 1: directly in the terminal or in the bash script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before executing the sith command (e.g. stretching), define the gaussian
function directly in the terminal or in the .sh script you want to run.

.. code-block:: bash

  # load your modules. For example,
  module load Gaussian/16.C.01-AVX2

  # Then, define the next function with your specific gaussian version.
  # In this case, we use g16 as an example.
  gaussian () {
    if [ -z "$2" ]
    then
      output="${1%.com}".log
    else
      output=$2
    fi

    g16 < $1 > $output
  }

  # source your sith command, in this case 'stretching'
  source $(sith stretching -path) -p 8 -i 1,37 -m molecule.pdb

^^^^^^^^^^^^^^^^^^^^^^^^^^^
Option 2: create a template
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The option 1 implies that every time you execute a sith command that requieres
gaussian, you have to run those lines that are probably always the same. To
avoid this redundancy, sith can be submitted directly as job script. For
that, you need to create the file :code:`$HOME/sw/load_modules.sh` with the
lines showed in the option 1, before sourcing sith. For example: 


.. code-block:: bash

  mkdir -p $HOME/sw
  cat << EOF > $HOME/sw/load_modules.sh
  # load your modules. For example,
  module load Gaussian/16.C.01-AVX2

  # Then, define the next function with your specific gaussian version.
  # In this case, we use g16 as an example.
  gaussian () {
    if [ -z "$2" ]
    then
      output="${1%.com}".log
    else
      output=$2
    fi

    g16 < $1 > $output
  }

Once this is done, you can use sith with the flag -c
(:code:`sith <function> -c ...`), for example:

.. code-block:: bash

  sith stretching -c -p 8 -i 1,37 -m molecule.pdb

Or you can even submit the sith command directly as a slurm job:




.. note::

  We designed sith to run in
  `slurm workload manager <https://slurm.schedmd.com/documentation.html>`_,
  which means you can submit your jobs directly with

  .. code-block:: bash

    sbatch -n 8 -J test $(sith stretching -path) -c -i 1,37 -m molecule.pdb

  To check the predefined options, run
  :code:`head -n 8 $(sith <function> -path)`. If you use any other queue
  manager, ignore this note.
