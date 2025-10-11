
=================
swap_atoms_in_com
=================

.. container:: bash-script-title

   :ref:`[script] <swap_atoms_in_com>` **sith/g09_stretching/swap_atoms_in_com.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Permutes the incides of two atoms in com file where the molecule is defined as
      a z-matrix.
      
        -a  <atom1> index of one of the atoms to be permuted. 1-based numbering.
        -b  <atom2> index of the other atom to be permuted. 1-based numbering.
        -f  <file.com> input gaussian file with the z-matrix already defined.
      
        -v  verbose.
        -h  prints this message.
      
