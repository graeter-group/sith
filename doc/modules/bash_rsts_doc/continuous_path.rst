
===============
continuous_path
===============

.. container:: bash-script-title

   :ref:`[script] <continuous_path>` **sith/g09_stretching/from_extreme/continuous_path.sh**

.. container:: bash-script-doc

   .. line-block::
      Creates the com files from the a set of xyz structures, and submit the jobs
      with :bashscript: (-S here
      refers to these subjobs)
      
        -c  Use this flag to run in a cluster. When -p is not defined, and you run in
            a slurm job manager system, the number of processors is equal to the
            number of cores asked in the submission of the job.
        -i  <index1,index2> indexes of the atoms used for constraining the distance
            of the intermedia structures when optimizing. If this flag is not used
            and a pdb is given through flag '-t', indexes 1 and 2 will correspond to
            the CH3 atoms in ACE and NME residues defined in the pdb if they exist.
        -l  <xc,base="bmk,6-31+g"> level of DFT theory.
        -n  <name> this script will collect all the *<name>*.xyz files sorted in
            alphabetic order as initial path. If this is a directory, this script
            will create a subdirectory there called 'conopt' and copy the
            *<name>*.xyz files to that subdirectory first.
        -p  <processors=1> number of processors per gaussian job. See description of
            flag -c.
        -S  <job_options=''> options for submitting a new job. This flag only makes
            sense in slurm cluster. Please, do not include a name (-J), nor the
            number of cores (-n, use -p for this). The input should be as in the next
            example: \"--partition=cpu --nice\".
        -t  <template.pdb> pdb file used to read the indexes if -i is not
            used.
      
        -v  verbose
        -h  prints this message.
