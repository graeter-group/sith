
==========
stretching
==========

.. container:: bash-script-title

   :ref:`[script] <stretching>` **sith/g09_stretching/stretching.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool obtains the stretched configurations of a molecule by increasing the
      distance between two atoms, constraining and optimizing at every step.
      
        -b  <number_of_breakages=1> The simulation will run until get this number of
            ruptures.
        -c  Use this flag to run in a cluster set with slurm. In that case, the
            of processors is equal to the number of cores asked in the submission of
            the job. When this flag is present, -p is ignored.
        -e  <extend_method=0> index of stretching method. To see the options, use
            'sith change_distance -h' to see the order.
            carbons of the capping groups
        -i  <index1,index2> indexes of the atoms to use for increasing the distance.
        -l  <xc,base=bmk,6-31+g> evel of DFT theory.
        -m  <molecule> molecule name. In this directory, a file called
            <molecule>-stretched00.pdb must exist.
        -p  <processors=1> number of processors per gaussian job.
        -r  restart stretching. In this case, this code must be executed from
            the molecule's directory.
        -s  <size[A]=0.2> Size of the step that increases the distances at each step.
      
        -v  verbose
        -h  prints this message.
      
