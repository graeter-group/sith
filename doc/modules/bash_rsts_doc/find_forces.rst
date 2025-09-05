
===========
find_forces
===========

.. container:: bash-script-title

   :ref:`[script] <find_forces>` **sith/g09_stretching/find_forces.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool computes the forces in all chk files and store them in a directory
      called forces.
      
        -c  run in cascade.
        -d  <dir=./> directory containging the chk files of the
            stretching-optimization process.
        -p  <pattern> pattern present in the chk files that will be used.
      
        -v  verbose.
        -h  prints this message.
      
      Note: it replaces the substring 'stretched' by 'forces' in the name.
      
