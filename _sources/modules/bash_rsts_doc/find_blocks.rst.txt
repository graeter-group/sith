
===========
find_blocks
===========

.. container:: bash-script-title

   :ref:`[script] <find_blocks>` **sith/utils/find_blocks.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code extracts the sections in a file starting and finishing with specific
      patterns without including the lines containing those patterns. Check the next
      options:
      
        -f  <file> file that shows 
        -s  <pattern> pattern that defines the beginning of the block. This line is
            not included in the block.
        -e  <pattern> pattern that defines the end of the block. This line is not
            included in the block.
        -i  use this flag if the start and the end are indexes
        -o  <output='output'> 'terminal' or the name of the output without extension.
            In the later case, the output will be stored in files called
            <output>_<n>.dat, where n is the number of appearence of the block in the
            file with 3-digits format (leading zeros); n starts in 001.
      
        -v  verbose of what's the code doing.
        -h  prints this message.
      
