
===========
restart_g09
===========



.. container:: bash-script-title

   :ref:`[script] <restart_g09>` **sith/g09_stretching/restart_g09.sh**

.. container:: bash-script-doc

   .. line-block::
      
      creates the com file to restart a gaussian09 job from an existing
      .chk file assuming that the com file of the process exists. I basically
      uses the same com file but replacing the geometry by the one in the .chk file.
      
        -c  run un a cluster.
        -d  <directory=./> directory where to find the .chk and .com files.
        -n  <name> name (without extension) of the .chk and .com files.
      
        -v  verbose.
        -h  prints this message.
      
