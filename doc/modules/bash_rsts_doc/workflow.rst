
========
workflow
========

.. container:: bash-script-title

   :ref:`[script] <workflow>` **sith/g09_stretching/workflow.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates all the stretched structures for a peptide and computes the
      needed quantities for sith. You can use this code to submit a Job in cluster or
      to execute it locally. Consider the next options:
      
        -b  <number of breakages=1> The simulation will run until getting this number
            of ruptures in the bonds.
        -c  run in cluster (see the documentation of the installation of sith
            -execute 'sith doc' in your terminal-)
        -i  <index1,index2> indexes of the atoms to use for increasing the distance.
            If these indices are not given and the molecule is an amino acid defined
            in a pdb, the CH3 atoms of the ACE and NME residues are chosen. 
        -l  <xc,base=bmk,6-31+g> evel of DFT theory.
        -m  <molecule> definition of the molecule (xyz, pdb, ...).
        -M  <method=0> Index of stretching method. To see the options, use
            'sith change_distance -h' to see the order.
        -n  <n_processors=1> number of processors per gaussian job.
        -r  restart. In this case, run from the directory of the pre-created
            peptide.
        -s  <size[A]=0.2> of the step that increases the distances.
        -S  <job_options=''> options for submitting a new job. This flag only makes
            sense in slurm cluster. Please, do not include a name and add the options
            as in the next example: "--partition=cpu --nice".
      
        -v  verbose.
        -h  prints this message.
      
