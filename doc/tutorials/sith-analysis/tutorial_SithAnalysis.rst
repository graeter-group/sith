.. _tutorial-SithAnalysis:

SITH Analysis Tutorial
======================

In this tutorial, you will learn how to analyse the outcome of your
stretching simulations with sith. You will obtain the distribution of
energies in a Glycine aminoacid with NME and ACE capping groups pulled
by the extremes (Figure 1). From the very beginning we assume that you have a
directory with `fchk` files in sith format with the QM information and you are
running the python scripts in this tutorial in the same directory. If you do
not have this files yet, download it from
`<https://github.com/Sucerquia/Gdata4tut>`_ or check the
:ref:`tutorial-stretching` if you want to generate it yourself.

.. code-block:: python
    :linenos:

    from sith.SITH import SITH

    sith = SITH('./forces')
    # The next line could also be jedi_analysis if your fchk files
    # contain the hessian matrix of the optimized configuration.
    sith.sith_analysis(); 


Visualize the Stretching of your Molecule
-----------------------------------------

.. note::

   `sith` has a visualization tool implemented, which requieres that you first
   install `VPython <https://vpython.org/>`_ and
   `VMol <https://github.com/Sucerquia/vmol>`_. **You can skip this part
   of the tutorial**.
   
   As an easy alternative to observe the stretching, you can use, for example,
   :code:`vmd *.xyz`.

First check the trajectory collected by :class:`sith.SITH.SITH`. It
comes from the fchk files you have in your directory ordered
alphabetically. In case of a disordered trajectory, rename the files in
increasing order.


.. code-block:: python
    :linenos:

    from sith.visualize.vmol import EnergiesVMol as Evmol
   
    v = Evmol(sith, absolute=True, deci=3,
              alignment=[1, 16, 13], height=300, default_bonds=True)


.. figure:: G-example_files/pulling-whole.gif
  :align: center

  **Figure 1.** Stretching path. The number shown in the box is the stretching
  state.

Notice that the dihedral angles of the hydrogens in the NME capping
group change drastically from the stretching **1** to the stretching
**2** (Figure 1), but we do not expect any distribution of energy in the degrees of
freedom of the hydrogen anyway because they are too light and not in the
connecting bonds (only side-chains). So, for the sake of the analysis
and to avoid numerical erros because of the large changes in the degrees
of freedom of the hydrogens, we can simply take them out of the anlysis.

.. code-block:: python
    :linenos:

    sith = SITH('./forces')
    sith.killer(killElements='H')
    sith.sith_analysis();

You can then visualize the evolution of the distribution of energies in
the interatomic distances (Figure 2).

.. code-block:: python
    :linenos:

    # The next two lines define the colormap used in this tutorial, 
    # you can use any of your preference instead, for example,
    # matplotlib.colormaps['Blues']
    import cmocean as cmo
    
    cmap = cmo.cm.algae 
    
    v = Evmol(sith, dofs=['bonds'], cmap=cmap, absolute=True, deci=4, portion=70,
              alignment=[1, 16, 13], labelsize=5)
    v.update_stretching(8)
    v.add_image_to_canvas(v.fig);

.. figure:: G-example_files/energy_in_bonds.png
  :align: center

  **Figure 2.** Energy distribution in distances at the 8th stretching step
  (right before the rupture of the molecule). The links to the hydrogen
  atoms do not appear because we toook them out of the analysis.

In the Figure 2, it is clear that the energy of the bonds at the 8th stretching
state are mostly stored in the C\ :math:`\alpha`-C bond of the ACE capping group,
followed by the C-N bond of the NME capping group. Of course, this is just a
visual inspection. In the next section, we run a more quantitative analysis of
this distribution of energies.

In the same way, you could view the distribution of energies only in the
angles or in the dihedrals. Or you could even include a comparison
considering angles and bonds (Figure 3).

.. code-block:: python
    :linenos:

    v = Evmol(sith, dofs=['angles', 'bonds'], cmap=cmap, absolute=True, deci=4, 
              alignment=[1, 16, 13], labelsize=5, default_bonds=True)
    v.update_stretching(8)
    v.add_image_to_canvas(v.fig);

.. figure:: G-example_files/energy_in_bonds_angles.png
  :align: center

  **Figure 3.** Energy distribution in distances and angles at the 8th
  stretching step (right before the rupture of the molecule). The links to the
  hydrogen atoms are gray because we toook them out of the analysis, but we
  added them as "default_bonds" for illustration proporsals.

Now, we can immediately observe that the energy is mostly stored in the
bonds, even when the energy loaded in some of the angles is higher that
then energy of some bonds.

.. note::

    The visualization of dihedrals might fail because of an internal error
    of vpython. Probably the error that you will see is

    .. parsed-literal::

        TypeError: unsupported operand type(s) for \*: 'float' and 'vpython.cyvector.vector'.



Quantitative Analysis
---------------------


Get to Know SITH Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the next lines of python commands, you could see the values stored
in the sith object. The name of the most important values are self
explanatory (just remember that dofs stands for Degrees of Freedom). For
further explanation of those values, check :class:`sith.SITH.SITH`.
For more clarity, here we discuss shortly how they help analyzing the
distribution of energies.

.. code-block:: python
    :linenos:

    import inspect
    
    [print(name) for name, value in inspect.getmembers(sith)
     if not callable(value) and not name.startswith("__")];


.. parsed-literal::

    all_dofs
    all_forces
    delta_q
    dim_indices
    dims
    dofs_energies
    energies_percertage
    n_structures
    reference
    removed_dofs
    structure_energies
    structures
    structures_scf_energies


In sith, all the quantities that depends on the stretching and the
degree of freedoms are numpy arrays where axis 0 and 1 run over the
stretching state and the degrees of freedom, respectively. The
stretching state is in increasing order and the index order of the
degrees of freedom is defined in :attr:`sith.SITH.SITH.dim_indices`.
Notice as well that :attr:`sith.SITH.SITH.dims` has 4 components that
correspond to the total number of Degrees of Freedom, the number of
distances, the number of angles and the number of dihedrals, in that
order.

SITH analysis
~~~~~~~~~~~~~

In this section we plot some of the important quantities in sith and
comment a bit on it.

.. note::

   We are using
   `PlottingTool <https://github.com/Sucerquia/PlotterTool/tree/master>`_
   and matplotlib, but feel free to use your desired plotting tool. The
   important lines are highlighted, the rest is just visualization
   preferences.

Internal Forces
^^^^^^^^^^^^^^^

In the next the python script, we plot the forces in the different
degrees of freedom (distances, angles and dihedrals). The important
attribute here is :attr:`sith.SITH.SITH.all_forces`. Notice that we plot
all the forces in all the stretching states (‘:’ as the selection of the
first component of all_forces -axis 0-), and we select the distances
first, then the angles, and finally, the dihedrals. Each plotted curve
correspond to the force at the stretching states in a given degree of
freedom labeled by the index of the DOF.

.. code-block:: python
    :linenos:
    :emphasize-lines: 10, 14, 18

    from PlottingTool import StandardPlotter
    import numpy as np
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(3, 1)
    sp = StandardPlotter(fig=fig, figheight=14,
                         ax=axes, ax_pref={'xticks': range(0, 9)})

    sp.plot_data(range(0, 9),
                 sith.all_forces[:, :sith.dims[1]].T,
                 ax=0, ms=2, lw=1,
                 data_label= np.arange(sith.dims[1]))
    sp.plot_data(range(0, 9),
                 sith.all_forces[:, sith.dims[1]: sith.dims[1] + sith.dims[2]].T,
                 data_label=np.arange(sith.dims[1], sith.dims[1] + sith.dims[2]),
                 ax=1, ms=2, lw=1)
    sp.plot_data(range(0, 9),
                 sith.all_forces[:, -sith.dims[3]:].T,
                 data_label= np.arange(sith.dims[0] - sith.dims[3], sith.dims[0]),
                 ax=2, ms=2, lw=1)
    
    sp.axis_setter(ax=0, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='F in Distance\n[Ha/Å]')
    sp.axis_setter(ax=1, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='F in Angle \n[Ha/rad]')
    sp.axis_setter(ax=2, legend=True, xlabel='Stretched Configuration',
                   ylabel='F in Dihedral\n[Ha/rad]')
    
    sp.spaces[0].set_axis(rows_cols=(3, 1), borders=[[0.16, 0.065], [0.99, 0.97]],
                          spaces=(0.1, 0.05));



.. figure:: G-example_files/G-example_23_0.png
    :align: center

    **Figure 4.** Forces per stretched structure for each degree of freedom
    separated by distances angles and dihedrals. The labels are the indices
    of the DOFS in sith.


Changes in the Degrees of Freedom
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other important attribute of sith is the change of the values in degrees
of freedom (:attr:`sith.SITH.SITH.delta_q`), which in this case is
defined by
:meth:`sith.energy_analysis.sith_analysis.SithAnalysis.get_sith_dq`.
This an important quantity (along with the forces) because of the
numerical integration, it is requiered to be small steps to conserve a
reasonable discretization of the energy integral. In this case, all the
changes in distances are lower than 0.1 Angstrom and the angles are
lower than 4 degrees.

.. code-block:: python
    :linenos:
    :emphasize-lines: 5, 9, 13

    fig, axes = plt.subplots(3, 1)
    sp = StandardPlotter(fig=fig, ax=axes, figheight=14)
    
    sp.plot_data(range(0, 9),
                 sith.delta_q[:, :sith.dims[1]].T,
                 data_label= np.arange(sith.dims[1]),
                 ax=0, ms=2, lw=1)
    sp.plot_data(range(0, 9),
                 sith.delta_q[:, sith.dims[1]: sith.dims[1] + sith.dims[2]].T * 180 / 3.14159,
                 data_label=np.arange(sith.dims[1], sith.dims[1] + sith.dims[2]),
                 ax=1, ms=2, lw=1)
    sp.plot_data(range(0, 9),
                 sith.delta_q[:, -sith.dims[3]:].T  * 180 / 3.14159,
                 data_label= np.arange(sith.dims[0] - sith.dims[3], sith.dims[0]),
                 ax=2, ms=2, lw=1)
    
    sp.axis_setter(ax=0, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='$\Delta$ Distances[Å]')
    sp.axis_setter(ax=1, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='$\Delta$Angles[\u00B0]')
    sp.axis_setter(ax=2, xminor=range(0, 9), legend=True, xlabel='Stretched Configuration',
                   ylabel='$\Delta$ Dihedrals[\u00B0]')
    
    sp.spaces[0].set_axis(rows_cols=(3, 1), borders=[[0.16, 0.1], [0.99, 0.99]],
                          spaces=(0.1, 0.05));



.. figure:: G-example_files/G-example_26_0.png
    :align: center

    **Figure 5.** :math:`\Delta q_i` per stretched structure for each degree of
    freedom separated by distances angles and dihedrals(:math:`q_i` is the
    value of the i-th degree of freedom). The labels are the indices of the
    DOFS in sith.

Values of the DOFs
^^^^^^^^^^^^^^^^^^

Instead of plotting the changes, we can plot the values of the degrees
of freedom at each stretching state, which are stored in
:attr:`sith.SITH.SITH.all_dofs`.

.. code-block:: python
    :linenos:
    :emphasize-lines: 5, 9, 13

    fig, axes = plt.subplots(3, 1)
    sp = StandardPlotter(fig=fig, ax=axes, figheight=14)
    
    sp.plot_data(range(0, 9),
                 sith.all_dofs[:, :sith.dims[1]].T,
                 data_label= np.arange(sith.dims[1]),
                 ax=0, ms=2, lw=1)
    sp.plot_data(range(0, 9),
                 sith.all_dofs[:, sith.dims[1]: sith.dims[1] + sith.dims[2]].T * 180 / 3.14159,
                 data_label=np.arange(sith.dims[1], sith.dims[1] + sith.dims[2]),
                 ax=1, ms=2, lw=1)
    sp.plot_data(range(0, 9),
                 sith.all_dofs[:, -sith.dims[3]:].T  * 180 / 3.14159,
                 data_label= np.arange(sith.dims[0] - sith.dims[3], sith.dims[0]),
                 ax=2, ms=2, lw=1)

    sp.axis_setter(ax=0, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='Distances[Å]')
    sp.axis_setter(ax=1, xticks=[], xminor=range(0, 9), legend=True,
                   ylabel='Angles[\u00B0]')
    sp.axis_setter(ax=2, xminor=range(0, 9), legend=True, xlabel='Stretched Configuration',
                   ylabel='Dihedrals[\u00B0]')
    
    sp.spaces[0].set_axis(rows_cols=(3, 1), borders=[[0.16, 0.1], [0.99, 0.99]],
                          spaces=(0.1, 0.05));

.. figure:: G-example_files/G-example_29_0.png
    :align: center

    **Figure 6.** :math:`q_i` per stretched structure for each degree of
    freedom separated by distances angles and dihedrals (:math:`q_i` is the
    value of the i-th degree of freedom). The labels are the indices of
    the DOFS in sith.

In the plot of the changes of the dihedrals obtained before, we did not
seem to have large changes, but the last plot seems to show large jumps.
However, those jumps are just an artifact of the visualization because
of the periodicity of the angles (181\ :math:`^\circ` =
-179\ :math:`^\circ`). We show below a way to obtain a better plot that
avoids this kind of jumps.

The Distribution of Energies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main point of sith is to compute the distribution of energies. This
quantity is stored in :attr:`sith.SITH.SITH.dofs_energies`.

.. code-block:: python
    :linenos:
    :emphasize-lines: 4-5, 7-9, 13-14

    fig, axes = plt.subplots(3, 1)
    sp = StandardPlotter(fig=fig, ax=axes, figheight=16)
    
    sp.plot_data(np.cumsum(sith.delta_q[:, :sith.dims[1]].T, axis=1),
                 sith.dofs_energies[:, :sith.dims[1]].T,
                 ax=0, ms=2, lw=1, data_label= np.arange(sith.dims[1]))
    sp.plot_data(np.cumsum(sith.delta_q[:,sith.dims[1]: sith.dims[1] + sith.dims[2]].T * 180 / 3.14159,
                           axis=1),
                 sith.dofs_energies[:, sith.dims[1]: sith.dims[1] + sith.dims[2]].T,
                 data_label=np.arange(sith.dims[1], sith.dims[1] + sith.dims[2]),
                 ax=1, ms=2, lw=1)
    
    sp.plot_data(np.cumsum(sith.delta_q[:, -sith.dims[3]:].T * 180 / 3.14159, axis=1),
                 sith.dofs_energies[:, -sith.dims[3]:].T,
                 data_label= np.arange(sith.dims[0] - sith.dims[3], sith.dims[0]),
                 ax=2, ms=2, lw=1)
    
    sp.axis_setter(ax=0,
                   xlabel='$\Delta$ Distance [Å]',
                   ylabel='$\Delta$ Energy [Ha]',
                   legend=True)
    
    sp.axis_setter(ax=1,
                   xlabel='$\Delta$ Angles[\u00B0]',
                   ylabel='$\Delta$ Energy [Ha]',
                   legend=True)
                   
    sp.axis_setter(ax=2,
                   xlabel='$\Delta$ Dihedrals[\u00B0]',
                   ylabel='$\Delta$ Energy [Ha]',
                   xticks=[-0.005, 0, 0.005, 0.01],
                   legend=True)
    
    sp.spaces[0].set_axis(rows_cols=(3, 1), borders=[[0.16, 0.065], [0.99, 0.97]],
                          spaces=(0.1, 0.1));



.. figure:: G-example_files/G-example_33_0.png
    :align: center

    **Figure 7.** Energy distribution per stretched structure for each degree of
    freedom separated by distances angles and dihedrals. The labels are the
    indices of the DOFS in sith.


Notice here, that the degree of freedom with index 5 (C-O disrance in
Glycie) does not store any energy, but it even shrinks a bit. Same for
the degree of freedom number 8 (C-O-C angle of the ACE capping group).
This is completely normal and comes from the relaxation of the molecule
when it is stretched, which implies that the change in a degree of
freedom is not enough reason to conclude that it stores any energy; it
should also have a force different to zero associated to it. In other
words, the external force alters the potential energy surface such that
the position of the minimum energy of this degree of freedom changes,
keeping the same value of energy.

Another important point to clarify is that the energy stored in the
dihedrals seem completely chaotic, and it happens because we are
observing only numerical zeros, which are actually random. Check again
the changes in the dihedrals and the forces in the plots before.

Finally, the most important result of these plots is to notice that all
the energy stored in the degrees of freedom increase instead of
decrease. This makes completely sense, because we start from a minimum
of energy from where the energy should not decrease but always increase.

Energy in One Specific DOF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Supposse that you want the energy storted only in the distance that
links the atoms 9 and 12 (1-based convention). First, you have to find
it in :attr:`sith.SITH.SITH.dim_indices`. Once you localize it, you
get the index that correspond to that degree of freedom (4 in this case)
and then use it for plotting the energy as in the following example:

.. code-block:: python
    :linenos:
    :emphasize-lines: 2, 3

    sp = StandardPlotter()
    sp.plot_data(np.cumsum(sith.delta_q[:, 4]),
                 sith.dofs_energies[:, 4],
                 data_label='9-12 distance', ms=2, lw=1);
    sp.axis_setter(xticks=[0, 0.05, 0.1, 0.15],xlabel='$\Delta$ Distance [Å]', legend=True,
                   ylabel='$\Delta$ Energy [Ha]')
    
    sp.spaces[0].set_axis(borders=[[0.145, 0.14], [0.99, 0.95]]);



.. figure:: G-example_files/G-example_36_0.png
    :align: center

    **Figure 8.** Energy stored in the degree of freedom with index 4, which
    corresponds to the distance between the atoms 9 and 12.


SithPlotter
-----------

We added a tool in sith for data visualization that is practical for
analysis, for example, if we are interested in a specific bond
recognizable by the name of the atoms involved (as shown below). This
tool is based on the existence of a pdb file that contains information
of the molecule like residue number and name of the atoms. In principle,
it was created for analysis of amino acids, but we are working to make
it more general to any kind of molecule. For now you might be able to
use the tools shown here thinking that amino is actually any residue in
your pdb.

.. code-block:: python
    :linenos:


    from sith.utils.sith_plots import SithPlotter
    
    pdb = './G.pdb'
    plotter_sith = SithPlotter(sith, pdb)

:class:`sith.utils.sith_plots` has the attributes `amino_name` and
`amino_info`. The first contains the name of the residies and the
second contains the name of the atoms per residue and the indices of
those atoms (1-based indices).

.. code-block:: python
    :linenos:

    print(plotter_sith.amino_name)
    print(plotter_sith.amino_info)


.. parsed-literal::

    {1: 'ACE ', 2: 'GLY ', 3: 'NME '}
    {1: {'CH3': 1, '1HH3': 2, '2HH3': 3, '3HH3': 4, 'C': 5, 'O': 6}, 2: {'N': 7, 'H': 8, 'CA': 9, 'HA1': 10, 'HA2': 11, 'C': 12, 'O': 13}, 3: {'N': 14, 'H': 15, 'CH3': 16, '1HH3': 17, '2HH3': 18, '3HH3': 19}}


In this way, if you want to select the index of the :math:`C_\alpha` and
N atoms of the Glycine (2nd residue), you only have to use the next
line:

.. code-block:: python
    :linenos:


    print(plotter_sith.amino_info[2]['CA'], plotter_sith.amino_info[2]['C']) 


.. parsed-literal::

    9 12


And if you want to know which is the index associated with the distance
that separates these two atoms, you can use

.. code-block:: python
    :linenos:

    plotter_sith.index_dof(np.array([9, 12, 0, 0]))


.. parsed-literal::

    4

Energy in DOF Defined by Atom Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this way, you can plot the energy in this specific degrees of freedom
as shown before. Or you can use the tool le_dof_amino.

.. code-block:: python
    :linenos:

    sp = StandardPlotter(ax_pref={'xlabel': r'C$\alpha$-C Length [Å]',
                                  'ylabel': r'E$_{C_\alpha-C}$ [Ha]',
                                  'xticks': [1.55, 1.6, 1.65, 1.7]})
    l, e = plotter_sith.le_dof_amino(('CA', 'C'), 2)
    sp.plot_data(l, e)
    sp.spaces[0].locate_ax(borders=[[0.14, 0.130], [0.99, 0.95]])


.. figure:: G-example_files/G-example_47_1.png
    :align: center

    **Figure 9.** Energy stored in the degree of freedom that corresponds to
    the distance between the C\ :math:`\alpha` and C atoms of the Glycine (2nd
    residue of our system).

Checking Harmonicity in one DOF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other methods like JEDI (T. Stauch and A. Dreuw (2016), Chem. Rev. 116) assume
a harmonic approximation in the distribution of
energies per degree of freedom. We can show with sith how close the
profile of an specific degree of freedom is to a harmonich potential. As
an example, we consider the the C\ :math:`\alpha`-C distance that we
just plotted. If the distribution of energies correspond to a harmonic
approximation in each degree of freedom, the next plot should correspond
to a straight line.

.. code-block:: python
    :linenos:

    sp = StandardPlotter(ax_pref={'xlabel': r'C$\alpha$-C Length [Å]',
                                  'ylabel': r'$\sqrt{E_{C_\alpha-C}} [\sqrt{Ha}]$',
                                  'xticks': [1.55, 1.6, 1.65, 1.7]})
    sp.plot_data(l, np.sqrt(e))
    sp.spaces[0].locate_ax(borders=[[0.14, 0.130], [0.99, 0.95]])


.. figure:: G-example_files/G-example_49_1.png
    :align: center

    **Figure 10.** Root squared of the energy stored in the degree of freedom
    that corresponds to the distance between the C\ :math:`\alpha` and C atoms
    of the Glycine (2nd residue of our system).

The Error in SITH
~~~~~~~~~~~~~~~~~

An interesting value to consider in sith is the error compared with the
DFT total change of energy. We know from the theory that the summation
of the change of energy in the DOFs should correspond with the total
change of energy of the whole molecule. In this case, then, we show that
our distribution of energies gives quite consistent results.

.. code-block:: python

    plotter_sith.plot_error();


.. figure:: G-example_files/G-example_51_1.png
    :align: center

    **Figure 11.** Different meassurements of the error in sith. The reference
    value is the change of total DFT energy.


Better Representation of Changes in Angles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the previous section, we showed that a normal plot might not be the
best idea to show values of angles, because of the periodicity. Then, we
added a method to :class:`SithPlotter` to have a better visualization.
In the next plot, you have the values if the angles

.. code-block:: python
    :linenos:

    sp = plotter_sith.plot_angles(step=2);



.. figure:: G-example_files/G-example_53_0.png
    :align: center

    **Figure 12.** Representation of angles and changes in polar projection.
