
===============
resubmit_failed
===============

.. container:: bash-script-title

   :ref:`[script] <resubmit_failed>` **sith/g09_stretching/from_extreme/resubmit_failed.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Takes the gaussian log files given by -l flag, and if the log file does not
      report a proper termination, a new job is resubmitted (creating a backup first)
      using the command given by the flag -e.
      
        -d  <path='./'> directory where the log and com files are.
        -e  <exec='ERROR: keyword single_g09 not recognized as part of pkgdeveloper. Use 'pkgdeveloper -h' to see the options you can use. -c -f '> execution command to be
            resubmitted. This command is very important to be given inside of " such
            that it is understood as only one value.
        -f  <frozen=''> line to freeze dofs. eg: '2 5 F'.
        -c  <comfile> input file used to run the previous trial.
        -l  <logfile> log file used to run the previous trial.
        -j  <jobname=comfile without extension> name of the new Job to be resubmited.
      
        -v  verbose
        -h  prints this message.
      
      Note: it assumes that ../frozen_dofs.dat exist.
      
      Note: In principle, this can be done easier by the chk, but then the input is
      not a Z-matrix anymore, then the output would not contain the info of the
      forces.
      
