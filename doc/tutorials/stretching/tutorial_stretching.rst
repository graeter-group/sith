.. _tutorial-stretching:

======================================
Energies Distribution in TriAlanine
======================================

.. attention::

    To use the tools of sith that are part of this tutorial, you have to
    import g09 as described in :ref:`load-modules`.

In this tutorial, we are going to extend a trialanine up to generate a rupture
in a bond. The extension is done by increasing the distance between two atoms,
fix the distance between those atoms and reoptimize. The process is repeated
with the optimized configuration. After obtaining each optimized configuration
we compute the forces in a set of internal degrees of freedom. With those
forces we compute the distribution of energies.

-------------------
Prepare your System
-------------------

Before you stretch your molecule using sith, you have to identify your molecule
and select the atoms that you want to use for the stretching. In this case we,
we are going to take a trialaline with ACE and NME capping groups. For the
stretching, we select the carbons of the methyl group in the capping ends. With
a visualization tool, you can check the index of those atoms (1-based indexing)
and see for yourself the indexes of the desired atoms are 1 and 16.
Alternatively, you can check these indexes by reading the pdb and looking for
the 'CH3' carbon atoms in the ACE and NME groups.

.. image:: selection.png
   :alt: Selection of the atoms to be pulled apart.
   :width: 640px
   :align: center

----------
Stretching
----------

Once you have your molecule and you know which atoms you want to pull apart,
you generate the pulling with the tool :code:`sith stretching`,

.. code-block:: bash

    sith stretching -c -i 1,16 -m G.pdb

.. note::

    For a whole description of stretching options, check **check** !!!

In 16 cores, this should take around 15 minutes to finish. The stretching
finishes after the cleavage of a bond. The broken molecule is stored in a
directory called rupture. In our case, the rupture happened between the methyl
group of the ACE capping and the rest of the molecule.

.. image:: rupture.png
   :alt: Broken molecule.
   :width: 640px
   :align: center

The outcome of the simulation is:

- A file called frozen_dofs.dat containing in the first row the atoms that
  used to pull the molecule with an F in front (in our case '1 16 F') and in
  the next line the atoms involved in the first rupture (in this case '1 5 F')
- A set of files called G-stretched<n>.<ext>, where <n> corresponds to the
  stretching step. <n>=0 means a free optimization without constraints. <ext>
  is:

  - chk: checkpoint file used by gaussian.
  - com: input file of the gaussian calculation. Notice that in the last line,
    the instruction of constraining the distance between the atoms is defined.
  - fchk: a human-readable chk file containing the final information of the
    optimization.
  - log: Gaussian log file.
  - xyz: coordinates of the optimized configurations.

- A directory called rupture containing the files of the molecule with the
  cleavage (the one showed above).






