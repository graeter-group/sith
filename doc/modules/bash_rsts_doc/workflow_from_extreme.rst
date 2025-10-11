
=====================
workflow_from_extreme
=====================

.. container:: bash-script-title

   :ref:`[script] <workflow_from_extreme>` **sith/g09_stretching/from_extreme/workflow_from_extreme.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates the files to do the sith analysis by optimizing a molecule,
      then takes the intermedia steps and prepares a constrained optimization of each
      one of them with :bashscript:
      after creating intermedias states that guarantee continuity in the degrees of
      freedom. This script submits in parallel each one of those optimizations with
      :bashscript: (-S here refers
      to these subjobs).
      
        -a  <alias> new name of the xyz file and subsequent files in directory
            'from_extreme'.
        -c  Use this flag to run in a cluster. When -p is not defined, and you run in
            a slurm job manager system, the number of processors is equal to the
            number of cores asked in the submission of the job.
        -i  <index1,index2> indexes of the atoms used for constraining the distance
            of the intermedia structures when optimizing. If this flag is not used
            but a pdb file is given (-t), indexes 1 and 2 will correspond to the CH3
            atoms in ACE and NME residues defined in the pdb if they exist.
        -l  <xc,base="bmk,6-31+g"> level of DFT theory.
        -m  <molecule> directory or coordinates file of configuration to be relaxed.
            For example, \"./AAA/\" a trialanine configuration ('last' after
            organizing all xyz files alphabetically all AAA*.xyz files in ./AAA/).
        -p  <processors=1> number of processors per gaussian job. See description of
            flag -c.
        -r  Use it to restart, in which case no directory will be created and the new
            <molecule>-optext.log is assumed to exist in the working directory.
        -S  <job_options=''> options for submitting a new job. This flag only makes
            sense in slurm cluster. Please, do not include a name (-J), nor the
            number of cores (-n, use -p for this). The input should be as in the next
            example: \"--partition=cpu --nice\".
        -t  <template.pdb> template pdb file to define indexes if -i is not used.
      
        -v  verbose.
        -h  prints this message.
      
      Note
      ----
      
        The outputs are stored in a directory called 'from_extreme'/
      
        This tool requires gaussian and sklearn.
