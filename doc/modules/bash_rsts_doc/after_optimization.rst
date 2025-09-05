
==================
after_optimization
==================

.. container:: bash-script-title

   :ref:`[script] <after_optimization>` **sith/g09_stretching/from_extreme/after_optimization.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Creates the com files from the xyz structures extracted from a gaussian log file
      and submit the corresponding jobs to compute the forces.
      
        -l  <log_file> optimization gaussian logfile.
        -n  <name> standard name. Usually pep name.
      
        -h  prints this message.
      
