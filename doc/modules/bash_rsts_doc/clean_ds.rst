
========
clean_ds
========

.. container:: bash-script-title

   :ref:`[script] <clean_ds>` **sith/g09_stretching/clean_ds.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code checks the aminoacids in the 'data runnig' directory and evaluates
      the state of the running (running, completed, failed). In case to be completed,
      it is moved to the dataset directory.
      
        -r  running directory. default ./
        -s  data set directory. default ../random3
        -u  user. defatult sucerqdl
      
        -v  verbose.
        -h  prints this message.
      
