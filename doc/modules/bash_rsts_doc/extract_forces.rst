==============
extract_forces
==============

.. container:: bash-script-title

   :ref:`[script] <extract_forces>` **sith/g09_stretching/extract_forces.sh**

.. container:: bash-script-doc

   .. line-block::

      Extract the forces and indexes of the DOFs from a log file (gaussian) and fchk
      file (or chk) when it exists. The output is a set of files called
      <pep>-forces<n_stretching>.fchk containing the information in fchk gaussian
      format.

        -f  <file.log> log file created by gaussian, chk file of this file should
            have the same name but different extension.
        -K  use this flag to avoid using chk file to extract forces.

        -v  verbose.
        -h   prints this message.
