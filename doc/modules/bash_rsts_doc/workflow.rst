
========
workflow
========

.. container:: bash-script-title

   :ref:`[script] <workflow>` **sith/g09_stretching/workflow.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates all the stretched structures for a peptide and computes the
      needed quantities for sith. You can use this code to submit a Job in cascade or
      to execute it locally. Consider the next options:
      
        -b  <number of breakages=1> The simulation will run until getting this number
            of ruptures in the bonds.
        -c  run in cascade. (modules are loaded)
        -d  <reference document=00-aminos.txt> file containing the existing peptides
            to avoid repetition. Mainly used for generation of random peptides.
        -e  <endo> or <exo> states for initial state of proline. Default 'random'.
        -m  <method=0> Index of stretching method. To see the options, use
            'sith change_distance -h' to see the order.
        -n  <options> Pepgen options (use \" for this)
        -p  <peptide> Chains of aminoacids to be evaluated. For example, "AAA"
            would analyse a trialanine peptide.
        -R  random pepeptide. Give the number of amino acids with this argument.
        -r  restart. In this case, run from the directory of the pre-created
            peptide.
        -s  <size[A]=0.2> of the step that increases the distances.
      
        -v  verbose.
        -h  prints this message.
      
