
=============
simplify_path
=============



.. container:: bash-script-title

   :ref:`[script] <simplify_path>` **sith/g09_stretching/simplify_path.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Takes all log files with certain pattern, extracts the last configuration,
      extract the dofs and removes the steps corresponding to oscillations or the
      system stacked in a configuration.
      
        -p  <pattern> pattern of the log files to be considered
      
        -v  verbose.
        -h  prints this message.
      
