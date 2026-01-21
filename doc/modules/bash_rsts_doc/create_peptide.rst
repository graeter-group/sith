
==============
create_peptide
==============



.. container:: bash-script-title

   :ref:`[script] <create_peptide>` **sith/g09_stretching/create_peptide.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates a peptide as a change of amino acids. It adds ACE and NME
      capping groups in the ends. The output is a pdb file with the name of the
      peptide saved in a directory with the same name.
      
        -d  <reference document=00-aminos.txt> file containing the existing peptides
            to avoid repetition. Mainly used for generation of random peptides.
        -e  <endo> or <exo> states for initial state of proline. Default 'random'.
        -n  <options> Pepgen options (use \" for this)
        -p  <peptide> Chains of aminoacids to be evaluated. For example, "AAA"
            would analyse a trialanine peptide.
        -R  random pepeptide. Give the number of amino acids with this argument.
      
        -v  verbose.
        -h  prints this message.
      
      Note
      ----
      
        This tools requires pepgen already installed.
      
