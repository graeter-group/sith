.. _tutorial-stretching:

======================================
Energies Distribution in Glycine
======================================

.. figure:: G-distribution.png
  :align: center

  **Figure 1.** Distribution of energies in a Glycine with ACE and NME capping
  groups pulled by an external force.

.. attention::

    To use the tools of sith that are part of this tutorial, you have to
    import g09 as described in :ref:`load-modules`.

In this tutorial, we are going to stretch a glycine amino acid with ACE and NME
capping groups from the carbon atoms in the ends of the peptide up to generate a
rupture in a bond (Figure 1). The extension is done by increasing the distance between two atoms,
fixing the distance between those atoms and reoptimizing. The process is repeated
with the optimized configuration. After obtaining each stretched configuration,
we compute the forces in a set of internal degrees of freedom. With those
forces we compute the distribution of energies.

-------------------
Prepare your System
-------------------

Before you stretch your molecule using sith, you have to identify your molecule
and select the atoms that you want to use for the stretching. In this case,
we are going to take a triglycine with ACE and NME capping groups. For the
stretching, we select the carbons of the methyl group in the capping ends
(Figure 1). With
a visualization tool, you can check the index of those atoms (1-based indexing)
and proof for yourself that the indexes of the desired atoms are 1 and 16.
Alternatively, you can check these indexes by reading the pdb file and looking for
the 'CH3' carbon atoms in the ACE and NME groups.

.. figure:: selection.png
  :align: center

  **Figure 1.** Illustration of the triglycine and the atoms used for the
  stretching.

----------
Stretching
----------

Once you have your molecule and you know which atoms you want to pull apart,
you generate the pulling with the tool
:bashscript:`sith.g09_stretching.stretching`

.. code-block:: bash

    sith stretching -c -i 1,16 -m G.pdb

In 16 cores, this should take around 15 minutes to finish; in one core it might
take around 2 hours. The stretching
finishes after the cleavage of a bond. The broken molecule is stored in a
directory called rupture. In our case, the rupture happened between the methyl
group of the ACE capping and the rest of the molecule (Figure 2).

.. figure:: rupture.png
  :align: center

  **Figure 2.** Molecule after cleavage.

The outcome of the simulation is:

- A file called frozen_dofs.dat containing in the first row the atoms
  used to pull the molecule with an F in front (in our case '1 16 F') and in
  the next line the atoms involved in the first rupture (in this case '1 5 F')
- A set of files called G-stretched<n>.<ext>, where <n> corresponds to the
  stretching step. <n>=0 means a free optimization without constraints. <ext>
  is:

  - chk: checkpoint file used by gaussian.
  - com: input file of the gaussian calculation. Notice the instruction of
    constraining the distance between the atoms is defined in the last line.
  - fchk: a human-readable chk file containing the final information of the
    optimization.
  - log: Gaussian log file.
  - xyz: coordinates of the optimized configurations.

- A directory called rupture containing the files of the molecule with the
  cleavage (the one showed in the Figure 2).

--------------
Extract forces
--------------

You just computed the optimized the stretched structures; you went through the
hardest part. Now you just need to compute the forces and set them in a proper
way that sith can read it. For this porpuse, use first
:bashscript:`sith.g09_stretching.find_forces`.

.. code-block:: bash

  sith find_forces -c -p "G-stretched" -v

This tool creates a directory called forces which contains a set of files
called G-forces<n>.<ext>, where <n> is the same as expained before and <ext> in
this case is com, log and chk (as described before). Finally, to have a proper
set of files that sith can read, we need to use 
:bashscript:`sith.g09_stretching.extract_forces`.


.. code-block:: bash

  sith extract_forces -c -d "forces"

This last step creates the set of fchk files which contain the information that
sith needs to compute the distribution of forces formated in according to
:mod:`sith.readers.g09_reader`. These fchk files are the outcome of the
gaussian :code:`formchk -3 G-force<n>.chk` and the log file.

-------------
SITH Analysis
-------------

You have now all the ingridients to cook your SITH analysis! Load the data and
compute the distribution of energies with your desired method.

.. code-block:: python

  from sith.SITH import SITH

  sith = SITH('./forces')
  sith.sith_analysis();

Then, the :class:`sith.SITH.SITH` object will contain all the necessary
information for further analysis. You could find thorough details of sith
analysis in the tutorial :ref:`tutorial-SithAnalysis`.
