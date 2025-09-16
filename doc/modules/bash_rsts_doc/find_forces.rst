
===========
find_forces
===========

.. container:: bash-script-title

   :ref:`[script] <find_forces>` **sith/g09_stretching/find_forces.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool computes the forces from a chk files that contains a given structure.
      
        -c  run in cascade.
        -f  <file> chk file.
        -p  <pattern> pattern present in the chk files that will be replaced with the
            word 'forces'.
      
        -v  verbose.
        -h  prints this message.
      
      Note: it replaces the substring 'stretched' by 'forces' in the name.
      
