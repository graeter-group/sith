.. SITH documentation master file, created by
   sphinx-quickstart on Fri Aug 26 12:18:30 2022.

.. mermaid::
   :align: center

   graph TD
   node1["sith"]
   click node1 "modules/sith.html" _self
     node1 --> node2["energy_analysis"]
     click node2 "modules/sith.energy_analysis.html" _self
     node1 --> node3["g09_stretching"]
     click node3 "modules/sith.g09_stretching.html" _self
       node3 --> node4["from_extreme"]
       click node4 "modules/sith.g09_stretching.from_extreme.html" _self
     node1 --> node5["readers"]
     click node5 "modules/sith.readers.html" _self
     node1 --> node6["utils"]
     click node6 "modules/sith.utils.html" _self
     node1 --> node7["visualize"]
     click node7 "modules/sith.visualize.html" _self

Contents
========

.. toctree::
   :maxdepth: 3
   :name: contents
   
   intro
   install
   workflow/workflow
   modules/sith
   tutorials/tutorials


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
