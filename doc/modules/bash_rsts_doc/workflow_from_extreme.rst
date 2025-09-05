
=====================
workflow_from_extreme
=====================

.. container:: bash-script-title

   :ref:`[script] <workflow_from_extreme>` **sith/g09_stretching/from_extreme/workflow_from_extreme.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates the files to do the sith analysis by optimizing a molecule
      that was just about to get a first rupture, then takes the intermedia steps and
      find the internal forces. Consider the next options:
      
        -c  run in cascade. (modules are loaded)
        -p  <peptide>. directory or xyzfile of last conf. Chains of aminoacids to
            be evaluated. For example, "./AAA/" would optimize the last
            stretched a trialanine peptide (where last means after organizing
            alphabetically).
        -r  Use if the job corresponds to a restart, in which case, no directory will
            be created and the new <peptide>-optext.com is assumed to exist.
      
        -v  verbose.
        -h  prints this message.
      
      Note: it is assumed that the file of the last configuration is named as:
      <amino acids-code>-<description><number of stretching>.xyz
      
