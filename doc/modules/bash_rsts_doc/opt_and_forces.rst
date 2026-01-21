
==============
opt_and_forces
==============



.. container:: bash-script-title

   :ref:`[script] <opt_and_forces>` **sith/g09_stretching/from_extreme/opt_and_forces.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code submits a gaussian job (typically an optimization) and uses the
      output to compute the forces.
      
        -c  Use this flag to run in a cluster. When -p is not defined, and you run in
            a slurm job manager system, the number of processors is equal to the
            number of cores asked in the submission of the job.
        -f  <com file> name of the gaussian input file without extension (.com). The
            output is has the same name, but .log extension.
        -F  use this flag fo AVIOD force calculation after optimization.
        -p  <processors=1> number of processors per gaussian job. See description of
            flag -c.
        -P  <pattern='conopt'> pattern to be replaced by 'forces' in the name of the
            output files of the optimization when submitting the forces.
        -S  <job_options=''> options for submitting a new job. This flag only makes
            sense in slurm cluster. Please, do not include a name and add the options
            as in the next example: "--partition=cpu --nice".
        -r  if the optimization did not converge, restart it from the last
            geometry. 
      
        -v  verbose.
        -h  prints this message.
      
