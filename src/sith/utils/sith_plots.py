from PlottingTool import StandardPlotter
from sith.utils.peptides import PepSetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from sith import SITH
from typing import Union, Tuple
from sith.utils.analysis import SithAnalysis
import cmocean as cmo
from matplotlib.colors import TwoSlopeNorm


def plot_averages_per_pos(ener_per_pos, mean_per_pos, stdr_per_pos, aminos,
                          ylabel, ylim):
    """
    Plot the mean value with error bars (std) per position (rows), per
    amino acids (x axis). It also shows the data points.

    Parameters
    ==========
    ener_per_pos: list
        set of energy data per position, per amino acid.
    mean_per_pos: list
        set of means per position, per amino acid.
    stdr_per_pos: list
        set of std per position, per amino acid.
    aminos: list
        list of amino acids in the data set.
    ylabel: str
        label of the y axis.
    ylim: list or tuple
        min and max value of the y axis.

    Return
    ======
    (StandardPlot) sith.plotters.StandardPlot object.

    Note
    ====
    The three first arguments can be obtained using analysis_energy_per_dof
    """
    indexes = [i for i in range(len(aminos))]

    fig, ax = plt.subplots(3, 1, figsize=(7, 7))
    sp = StandardPlotter(fig=fig, ax=ax,
                         ax_pref={'xticks': indexes,
                                  'xticklabels': [None] * len(indexes),
                                  'xlim': [indexes[0] - 1, indexes[-1] + 1],
                                  'grid': True,
                                  'color_grid': [0.95, 0.95, 0.95],
                                  'ylim': ylim})

    for i in [0, 1, 2]:
        for j, _ in enumerate(ener_per_pos[i]):
            y = ener_per_pos[i][aminos[j]]
            if y is None:
                continue
            x = [j + (np.random.rand() - 0.5) * 0.2 for _ in y]
            sp.plot_data(x,
                         y, ax=i, pstyle='o', markersize=2,
                         color_plot=[0.7, 0.7, 0.7], fillstyle='none')
            sp.ax[i].errorbar(j, mean_per_pos[i][aminos[j]],
                              yerr=stdr_per_pos[i][aminos[j]],
                              fmt='o', color='C0',
                              markersize=2, lw=0.7)

        props = dict(boxstyle='round', facecolor='white', lw=0.5)
        factor = 0.92
        sp.ax[i].text(-0.5, ylim[-1] * factor + (1 - factor) * ylim[0],
                      f"Position {i+1}", fontsize=8, va='top',
                      weight='bold', color=[0.4, 0.4, 0.4], bbox=props)
        sp.ax[i].spines[['right', 'top']].set_visible(False)

    sp.axis_setter(ax=2, xlabel="Amino Acid", xticklabels=aminos)
    sp.axis_setter(ax=1, ylabel=ylabel)

    sp.spaces[0].set_axis(borders=[[0.105, 0.08], [0.99, 0.97]],
                          spaces=(0.1, 0.07),
                          rows_cols=(3, 1))

    return sp


def plot_matrix(matrix, labels, n_per_ele, cbar_label):
    """
    Plot a matrix as a heatmap and show a value per combination.

    Parameters
    ==========
    matrix: numpy.array
        squared matrix as a shape of matrix of dimension 2.
    labels: list
        label of the components of the matrix.
    n_per_ele: numpy.array
        value to be displayed on the box of the matrix. For example, the number
        of samples.
    cbar_label: str
        label of the color bar.

    Return
    ======
    (StandardPlot) sith.plotters.StandardPlot object.
    """
    sp = StandardPlotter(figwidth=8.9, figheight=8.9 / 1.3,
                         ax_pref={'xticks': np.arange(matrix.shape[1]),
                                  'xticklabels': labels,
                                  'yticks': np.arange(matrix.shape[0]),
                                  'yticklabels': labels})

    # Plot the heatmap
    cmap = cmo.cm.algae
    im = sp.ax[0].imshow(matrix, cmap=cmap)

    # Setting matrix plot
    sp.ax[0].set_yticks(np.arange(matrix.shape[0]), labels=labels)
    sp.ax[0].tick_params(top=True, bottom=False, labeltop=True,
                         labelbottom=False)
    sp.ax[0].set_xticks(np.arange(matrix.shape[1] + 1) - 0.5, minor=True)
    sp.ax[0].set_yticks(np.arange(matrix.shape[0] + 1) - 0.5, minor=True)
    sp.ax[0].grid(which="minor", color=[0.4, 0.4, 0.4], linestyle='--',
                  linewidth=0.3)
    sp.ax[0].tick_params(which="minor", bottom=False, left=False)
    sp.ax[0].set_xlabel('Pos 1 or 3', labelpad=-190)
    sp.ax[0].set_ylabel('Pos 2', labelpad=0)

    # Create colorbar
    sp.add_axes()
    # remove next line and add reset in the setter

    plt.colorbar(im, cax=sp.ax[1])
    sp.axis_setter(ax=sp.ax[1], ylabel=cbar_label, reset=True)
    sp.ax[1].yaxis.get_offset_text().set_x(1.55)

    # n_numbers
    for i in range(len(labels)):
        for j in range(len(labels)):
            color = "black"
            if matrix[j][i] > (np.nanmax(matrix.flatten())
                               + np.nanmin(matrix.flatten())) / 2:
                color = "white"
            sp.ax[0].text(i, j, int(n_per_ele.T[i][j]), va='center',
                          ha='center', fontsize=5, color=color)

    # ==== Location of axis
    # Matrix:
    side = 0.87
    left = -0.005
    bottom = -0.04
    sp.spaces[0].locate_ax([[left + 0.0, bottom + (1 - side) / 2],
                            [left + side, bottom + side]],
                           ax=sp.ax[0])
    # ColorBar
    sp.spaces[0].locate_ax([[0.76, bottom + (1 - side) / 2], [0.82, 0.912]],
                           ax=sp.ax[1])

    return sp


def plot_matrix2(matrix, labels, cbar_label):
    """
    Plot a matrix as a heatmap. The colormap is plt.get_cmap('coolwarm') and is
    divergent, with 0.05 in the middle.

    Parameters
    ==========
    matrix: numpy.array
        squared matrix as a shape of matrix of dimension 2.
    labels: list
        label of the components of the matrix.
    cbar_label: str
        label of the color bar.

    Return
    ======
    (StandardPlot) sith.plotters.StandardPlot object.

    Note
    ====
    Usually applied to show p-values.
    """
    sp = StandardPlotter(figwidth=8.9, figheight=8.9 / 1.3,
                         ax_pref={'xticks': np.arange(matrix.shape[1]),
                                  'xticklabels': labels,
                                  'yticks': np.arange(matrix.shape[0]),
                                  'yticklabels': labels})

    # Plot the heatmap
    cmap = plt.get_cmap('coolwarm')
    norm = TwoSlopeNorm(vmin=0, vcenter=0.05, vmax=1)
    im = sp.ax[0].imshow(matrix, cmap=cmap, norm=norm)

    # Setting matrix plot
    sp.ax[0].set_yticks(np.arange(matrix.shape[0]), labels=labels)
    sp.ax[0].tick_params(top=True, bottom=False, labeltop=True,
                         labelbottom=False)
    sp.ax[0].set_xticks(np.arange(matrix.shape[1] + 1) - 0.5, minor=True)
    sp.ax[0].set_yticks(np.arange(matrix.shape[0] + 1) - 0.5, minor=True)
    sp.ax[0].grid(which="minor", color=[0.4, 0.4, 0.4], linestyle='--',
                  linewidth=0.3)
    sp.ax[0].tick_params(which="minor", bottom=False, left=False)

    # Create colorbar
    sp.add_axes()

    sp.ax[-1].preferences = sp.ax_pref_bck
    plt.colorbar(im, cax=sp.ax[1])
    sp.axis_setter(ax=sp.ax[1], ylabel=cbar_label, reset=True)
    sp.ax[1].yaxis.get_offset_text().set_x(1.55)

    # ==== Location of axis
    # Matrix:
    side = 0.87
    left = -0.005
    bottom = -0.04
    sp.spaces[0].locate_ax([[left + 0.0, bottom + (1 - side) / 2],
                            [left + side, bottom + side]],
                           ax=sp.ax[0])
    # ColorBar
    sp.spaces[0].locate_ax([[0.76, bottom + (1 - side) / 2], [0.82, 0.912]],
                           ax=sp.ax[1])

    return sp


def plot_matrix3(matrix, labels, cbar_label):
    """
    Plot a matrix as a heatmap and show a value per combination.

    Parameters
    ==========
    matrix: numpy.array
        squared matrix as a shape of matrix of dimension 2.
    labels: list
        label of the components of the matrix.
    n_per_ele: numpy.array
        value to be displayed on the box of the matrix. For example, the number
        of samples.
    cbar_label: str
        label of the color bar.

    Return
    ======
    (StandardPlot) sith.plotters.StandardPlot object.
    """
    sp = StandardPlotter(figwidth=8.9, figheight=8.9 / 1.3,
                         ax_pref={'xticks': np.arange(matrix.shape[1]),
                                  'xticklabels': labels,
                                  'yticks': np.arange(matrix.shape[0]),
                                  'yticklabels': labels})

    # Plot the heatmap
    cmap = cmo.cm.algae
    im = sp.ax[0].imshow(matrix, cmap=cmap)
    # Setting matrix plot
    sp.ax[0].set_yticks(np.arange(matrix.shape[0]), labels=labels)
    sp.ax[0].tick_params(top=True, bottom=False, labeltop=True,
                         labelbottom=False)
    sp.ax[0].set_xticks(np.arange(matrix.shape[1] + 1) - 0.5, minor=True)
    sp.ax[0].set_yticks(np.arange(matrix.shape[0] + 1) - 0.5, minor=True)
    sp.ax[0].grid(which="minor", color=[0.4, 0.4, 0.4], linestyle='--',
                  linewidth=0.3)
    sp.ax[0].tick_params(which="minor", bottom=False, left=False)

    # Create colorbar
    sp.add_axes()

    plt.colorbar(im, cax=sp.ax[1])
    sp.axis_setter(ax=sp.ax[1], ylabel=cbar_label, reset=True)
    sp.ax[1].yaxis.get_offset_text().set_x(1.55)

    # ==== Location of axis
    # Matrix:
    side = 0.87
    left = -0.005
    bottom = -0.04
    sp.spaces[0].locate_ax([[left + 0.0, bottom + (1 - side) / 2],
                            [left + side, bottom + side]],
                           ax=sp.ax[0])
    # ColorBar
    sp.spaces[0].locate_ax([[0.76, bottom + (1 - side) / 2], [0.82, 0.912]],
                           ax=sp.ax[1])

    return sp


class SithPlotter(PepSetter, SithAnalysis):
    """
    Object that plots the main graphs to analyze sith outcomes.
    """
    def __init__(self, sith: SITH, pdb_template: str):
        """
        Parameters
        ==========
        sith:
            sith object containing all information about the sith analyzis.
        pdb_remplate:
            path to .pdb file that has the peptide information.
        """
        self.sith = sith
        PepSetter.__init__(self, pdb_template)
        SithAnalysis.__init__(self, self.sith, self)

    def plot_energies_in_DOFs(self, steps: list = None,
                              jump_stretching: int = 1,
                              **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot of distribution of energies in all degrees of freedom and in each
        kind. Namely, distances, angles, dihedrals. Then, it creates a 2x2
        plot.

        Parameters
        ==========
        steps: list. Default=[1, 1, 1, 1]
            size of steps separating the labels of the degrees of freedom.
        **kwargs:
            SithPlotter.plot_data arguments.

        Return
        ======
        plt.figure.Figure, plt.Axes
        """
        if steps is None:
            steps = [1, 1, 1, 1]
        fig, axes = plt.subplots(2, 2)
        sp = StandardPlotter(fig=fig, ax=axes, ax_pref={'ylabel': '',
                                                        'l_ticks': 2})
        plots_space = sp.add_space(borders=[[0, 0], [0.825, 1]])
        plots_space.set_axis(rows_cols=(2, 2), borders=[[0.13, 0.13],
                                                        [1, 0.95]],
                             spaces=(0.07, 0.185))
        energies_per_DOF = self.sith.dofs_energies
        dims = self.sith.dims

        emin = min(energies_per_DOF.flatten())
        emax = max(energies_per_DOF.flatten())

        # Add separation of dofs
        sp.plot_data([dims[1] + 0.5, dims[1] + 0.5], [emin, emax],
                     pstyle='--', color_plot='gray', ax=0)
        sp.plot_data([dims[1] + dims[2] + 0.5, dims[1] + dims[2] + 0.5],
                     [emin, emax], pstyle='--', color_plot='gray', ax=0)

        # Plot energies
        self.plot_sith(np.arange(1, dims[0] + 1), energies_per_DOF,
                       'All DOFs', ax=0, sp=sp, cbar=False, step=steps[0],
                       show_amino_legends=True,
                       jump_stretching=jump_stretching,
                       **kwargs)
        self.plot_sith(np.arange(1, dims[1] + 1),
                       energies_per_DOF[:, :dims[1]],
                       'Distances', ax=1, cbar=False, step=steps[1], sp=sp,
                       jump_stretching=jump_stretching, **kwargs)
        self.plot_sith(np.arange(dims[1] + 1, dims[1] + dims[2] + 1),
                       energies_per_DOF[:, dims[1]:dims[1] + dims[2]],
                       'Angles', ax=2, cbar=False, step=steps[2], sp=sp,
                       jump_stretching=jump_stretching, **kwargs)
        self.plot_sith(np.arange(dims[1] + dims[2] + 1, dims[0] + 1),
                       energies_per_DOF[:, dims[1] + dims[2]:],
                       'Dihedrals', ax=3, cbar=True,
                       step=steps[3], sp=sp, jump_stretching=jump_stretching,
                       **kwargs)

        # Remove Unnecessary labels
        for i in [1, 3]:
            sp.axis_setter(ax=i, ylabel='')

        return sp

    def plot_sith(self, dofs: Union[list, tuple, np.ndarray] = None,
                  e_dofs: Union[list, tuple, np.ndarray] = None,
                  xlabel: str = '', ax: Union[plt.Axes, int] = 0,
                  sp: StandardPlotter = None,
                  cmap: mpl.colors.Colormap = None,
                  cbar: bool = True, step: int = 1, pstyle: str = '-o',
                  ylabel: str = r'$\Delta$E$_{\rm{\bf i}}$[' + f'Ha]',
                  jump_stretching: int = 1,
                  show_amino_legends: bool = False, ax_pref: dict = {},
                  sp_pref={},
                  pad_cbar=10,
                  **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        This function plots the energies per degrees of freedom from
        SithPlotter.sith.energies

        Parameters
        ==========
        e_dofs: array
            labels of the degrees of freedom.
        dofs: array
            energies per degree of freedom. Usually a matrix where each
            component contains the energies for each DOF for each deformed
            config:
            dof\\ deformed     0 1  2  3 ...
            0             [[              ]]
            1             [[              ]]
            2             [[              ]]
            .
            .
            .
        xlabel: str
            label of the xlabel indicating the represented DOFS.
        sp: StandardPlotter
            plotter object. if not given. It creates a new object with one
            graph.
        cmap: plt.color.Colormap. Default: cmocean.cm.algae or 'vidris'.
            Color map for the deformations.
        cbar: bool. Default=False
            True to show the color bar.
        step: int. Default 1
            size of steps separating the labels of the degrees of freedom.
        pstryle: str. Default='-o'
            style of the lines of energies
        ylabel: str. Default={r'$\Delta$E$_{\\rm{\bf i}}$' + f'Ha'}
            label in the y axis.
        jump_stretching: int. Default=1
            jumps from one stretching to the other starting from the optimized
            one.
        show_amino_legends: bool. Default=False
            True to show the name of the aminoacids painting the background.
        ax_pref: dict
            axes preferences. See StandardPlotter.

        Return
        ======
        (plt.Figure, plt.Axes) plotting objects used to create the figure.
        """
        if dofs is None:
            dofs = np.arange(1, self.sith.dims[0] + 1)
        if e_dofs is None:
            e_dofs = self.sith.dofs_energies
        e_dofs = e_dofs[::jump_stretching]

        # Setup default
        if cmap is None:
            try:
                import cmocean as cmo
                cmap = cmo.cm.algae
            except ImportError:
                cmap = mpl.get_cmap['viridis']

        if sp is None:
            sp = StandardPlotter(**sp_pref)

        if isinstance(ax, int):
            ax = sp.ax[ax]
        else:
            raise ValueError("\"ax\" must be an intiger in plot_sith method")

        fig = sp.fig

        # Color bar
        boundaries = np.arange(0, len(e_dofs) + 1)
        normalize = mpl.colors.BoundaryNorm(boundaries - 0.5, cmap.N)
        if cbar:
            ax_bar = sp.add_axes()
            space_bar = sp.add_space(borders=[[0.845, 0], [1, 1]],
                                     axes=ax_bar)
            cbar = sp.fig.colorbar(mpl.cm.ScalarMappable(norm=normalize,
                                                         cmap=cmap),
                                   cax=ax_bar,
                                   orientation='vertical')
            cbar.set_ticks(boundaries[:-1],
                           labels=boundaries[:-1] * jump_stretching)
            cbar.set_label(label="Stretched",
                           fontsize=mpl.rcParams['font.size']
                           * sp.ax_pref['labels_scale'],
                           rotation=90,
                           labelpad=pad_cbar)
            cbar.ax.tick_params(length=0)
            space_bar.locate_ax(borders=[[0, 0.13], [0.2, 0.95]])

        sp.axis_setter(ax=ax, xlabel=xlabel, ylabel=ylabel,
                       xticks=np.arange(dofs[0], dofs[-1] + 1, step),
                       **ax_pref)

        [sp.plot_data(dofs, e_dofs[i], ax=ax, pstyle=pstyle,
                      color_plot=cmap(normalize(i)),
                      **kwargs
                      ) for i in range(len(e_dofs))]

        colors = self.add_color_per_amino(ax)
        if show_amino_legends:
            ax.legend(handles=list(colors.values()), loc='upper right',
                      fontsize=mpl.rcParams['font.size']
                      * sp.ax_pref['ticks_scale'])
        ax.set_xlim([dofs[0] - 0.5, dofs[-1] + 0.5])

        return sp

    def add_color_per_amino(self, ax: plt.Axes) -> dict:
        """
        Add an colored rectangle in the background for every DOF belonging to
        an aminoacid.

        Paramenters
        ===========
        ax: plt.Axes
            Axes of the graphics to add the colors

        Return
        ======
        (dict) colors patches per amino acid labeled by indices.
        """
        dofs_classified = self._dof_classificator()
        cmap = plt.get_cmap('tab10_r')
        boundaries = np.arange(1, 11, 1)
        normalize = mpl.colors.BoundaryNorm(boundaries - 0.5, cmap.N)
        patches = {}
        for i in dofs_classified.keys():
            blocks = self.create_blocks(dofs_classified, i)
            for region in blocks:
                ax.axvspan(region[0] + 0.5,
                           region[1] + 1.5,
                           color=cmap(normalize(i)),
                           alpha=0.1)
                patch = mpatches.Patch(facecolor=cmap(normalize(i)),
                                       label=f'{i}-{self.amino_name[i]}',
                                       alpha=0.1,
                                       edgecolor="black", linewidth=1)
                patches[i] = patch
        return patches

    def _dof_classificator(self):
        """"
        classify the degrees of freedom according to the aminoacid del belong.

        Return
        ======
        (dict) The keys are the index of the amino acid, the values are the
        list of DOFs belonging to them.
        """
        atoms_per_aminoacids = self.atom_indexes
        dofs_indexes = self.sith.dim_indices
        dof_per_amino = {}
        for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
            dof_per_amino[i] = np.array([], dtype=int)

        for i in range(len(dofs_indexes)):
            for j in atoms_per_aminoacids.keys():
                if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]
                           + [0]).all():
                    dof_per_amino[j] = np.append(dof_per_amino[j], i)
                    break
        return dof_per_amino

    def create_blocks(self, classified, key):
        blocks = [[classified[key][0], classified[key][0]]]
        for i in classified[key][1:]:
            if i == blocks[-1][-1] + 1:
                blocks[-1][-1] = i
            else:
                blocks.append([i, i])
        return blocks

    def plot_angles(self,
                    cmap: mpl.colors.Colormap = None,
                    step: int = 1,
                    sp_pref: dict = {}) -> Tuple[plt.Figure,
                                                 plt.Axes]:
        """
        Plot values of angles and changes during the deformations.

        Parameters
        ==========
        cmap: Colormap
            colormap to the increasing changes.
        step: int
            steps between radius ticks.
        sp_pref: dict
            dictionary with kwargs for StandardPlotter.

        Return
        ======
        (StandardPlotter) Standard plotter.
        """
        # Double column format
        if 'figwidth' not in sp_pref:
            sp_pref['figwidth'] = 18.3

        if 'figheight' not in sp_pref:
            sp_pref['figheight'] = 15

        sith = self.sith
        distances = sith.dims[1]
        n_angles = sith.dims[2] + sith.dims[3]
        n_deformed = sith.n_structures
        rs = [np.arange(1, n_deformed + 1) for _ in range(n_angles)]

        fig, axes = plt.subplots(2, 2)
        sp = StandardPlotter(fig=fig, ax=axes, **sp_pref)
        scale = 0.03
        # polar plots
        sp.set_polar(ax=1, r_ticks=rs[0][::step],
                     r_lims=[-rs[0][-1] * scale, rs[0][-1] * (1 + scale)])
        sp.set_polar(ax=3, r_ticks=rs[0][::step],
                     r_lims=[-rs[0][-1] * scale, rs[0][-1] * (1 + scale)])
        sp.plot_data(sith.all_dofs[:, distances:].T,
                     rs,
                     pstyle='-', ax=1)
        sp.plot_data(sith.delta_q[:, distances:].T,
                     rs,
                     pstyle='-', ax=3)
        sp.spaces[0].set_axis(rows_cols=(2, 2),
                              borders=[[0.1, 0.08], [0.97, 0.97]],
                              spaces=(0.05, 0.1), axes=sp.ax)

        # Plot values in cartesian
        if cmap is None:
            try:
                import cmocean as cmo
                cmap = cmo.cm.algae
            except ImportError:
                cmap = mpl.colormaps['viridis']
        deformations = np.arange(1, n_deformed + 1, 1)
        normalize = mpl.colors.BoundaryNorm(deformations + 0.5, cmap.N)
        colors = [cmap(normalize(i + 0.5))[:3] for i in deformations]

        # Set axes
        sp.axis_setter(ax=0, xlabel='Angle index', ylabel='Value[rad]')
        sp.axis_setter(ax=2, xlabel='Angle index', ylabel='Changes[rad]')

        # Add limits at pi and -pi
        [sp.plot_data([distances + 1, distances + n_angles],
                      [[np.pi, np.pi], [-np.pi, -np.pi]],
                      pstyle='--', color_plot='gray', ax=i)
         for i in [0, 2]]

        sp.plot_data(np.arange(distances + 1, distances + n_angles + 1),
                     sith.all_dofs[:, distances:], pstyle='s',
                     markersize=3, color_plot=colors)
        sp.plot_data(np.arange(distances + 1, distances + n_angles + 1),
                     sith.delta_q[:, distances:], pstyle='s',
                     markersize=3, ax=2, color_plot=colors)

        # Fix location
        ho = 0.6
        # adjust cartesians
        sp.add_space(borders=[[0, 0], [ho, 1]],
                     axes=sp.ax[[0, 2]])
        sp.spaces[-1].set_axis(rows_cols=(2, 1),
                               borders=[[0.125, 0.075], [0.99, 0.99]],
                               spaces=(0.05, 0.1))
        # adjust polars
        sp.add_space(borders=[[ho, 0], [1, 1]],
                     axes=sp.ax[[1, 3]])
        sp.spaces[-1].set_axis(rows_cols=(2, 1),
                               borders=[[0.125, 0.05], [0.9, 0.95]],
                               spaces=(0.05, 0.1))

        return sp

    def plot_error(self,
                   classical: Union[list, tuple, np.ndarray] = None,
                   sp_pref={}) -> Tuple[np.ndarray, np.ndarray, list]:
        """
        Plot the error between the expected value (DFT) and the computed using
        sith.

        Parameters
        ==========
        classical: array-like
            set of classical energies of each deformation computed with
            amber99.
        sp_pref:


        Return
        ======
        (np.ndarray, np.ndarray, list) energies, distances and sp.axes
        """
        if 'figheight' not in sp_pref:
            sp_pref['figheight'] = 10
        # Check if ACE is the last or the first residue
        first_cap = self.amino_name[1]
        if first_cap == 'ACE':
            first_atom = 'N'
            last_atom = 'C'
        else:
            first_atom = 'C'
            last_atom = 'N'

        # Find the index of the closest atoms to the cap residues.
        last_amino = list(self.amino_info.keys())[-2]
        index1 = self.amino_info[2][first_atom]
        index2 = self.amino_info[last_amino][last_atom]

        # Find distances
        distances = []
        for defo in self.sith.structures:
            distances.append(defo.atoms.get_distance(index1 - 1, index2 - 1))
        distances = (np.array(distances) - distances[0])

        # get energies dE-sith, dE-DFT, dE-error, dE-errorpercent
        e = self.sith.energies_error()
        e = (np.sum(self.sith.dofs_energies, axis=1),
             self.sith.structures_scf_energies,
             e[0],
             e[1])

        # ==== Plot ====
        # set up
        fig, axes = plt.subplots(3, 1, figsize=(5, 13))
        sp = StandardPlotter(fig=fig, ax=axes, **sp_pref)
        ticks = np.round(np.linspace(0, distances[-1], 6), decimals=2)
        ws = (ticks[1] - ticks[0]) / 5
        sp.spaces[0].set_axis(rows_cols=(3, 1), spaces=(1, 0.03),
                              borders=[[0.2, 0.1], [0.98, 0.99]])
        sp.axis_setter(ax=0, ylabel='$\Delta$E [Ha]',
                       xticks=[], xminor=ticks,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        sp.axis_setter(ax=1,
                       ylabel='$\Delta$E$_{SITH}$ - $\Delta$E$_{DFT}$ [Ha]',
                       xticks=[], xminor=ticks,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        sp.axis_setter(ax=2, ylabel='Error [%]',
                       xlabel='$\Delta$ End-to-end distance [Å]',
                       xticks=ticks, xminor=ticks - 0.0001,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        # plot axis 1
        sp.plot_data(distances, e[1], ax=0, data_label='DFT', pstyle='-')
        sp.plot_data(distances, e[0], ax=0, data_label='SITH', pstyle='o',
                     fillstyle='none', ms=5)
        if classical is not None:
            sp.plot_data(distances, classical, ax=0, data_label='amber99',
                         pstyle='*-', fraclw=10)
        sp.ax[0].legend()

        # plot axis 2
        sp.plot_data(distances, e[2].T, ax=1, pstyle='o-', color_plot='C2')

        # plot axis 3
        sp.plot_data(distances, e[3].T, ax=2, pstyle='o-', color_plot='C2')

        sp.spaces[0].set_axis(rows_cols=(3, 1),
                              borders=[[0.15, 0.11], [0.99, 0.99]],
                              spaces=(0.05, 0.06))

        return e, distances, sp

    def plot_hessian(self, ax=None, deci=2, orientation='vertical', cbar=True,
                     ticks=15):
        """
        Function that plots the a matrix using a divergent colormap to separate
        the negative from the positive values.

        Parameters
        ==========
        hessian: NxN numpy.array
            matrix to be ploted
        ax: plt.Axes
            Axis to add the plot. Default: None, in this case, the function
            creates a new Axis.
        deci: int
            number of decimals in the colorbar.
        orientation: str
            orientation of the colorbar. Default: 'vertical'.
        cbar: Bool
            True to show the colorbar. Default: True
        ticks: float
            ticks size.

        Return
        ======
        PathCollection
        """
        ref = self.sith.reference
        hessian = self.sith.structures[ref].hessian

        if orientation[0] == 'v':
            pad = 0.02
            shrink = 1
            rotation = 0
        else:
            pad = 0.15
            shrink = 0.9
            rotation = 90

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(10, 10))
        if orientation[0] == 'v':
            pad = 0.02
            shrink = 0.85
            rotation = 0
        else:
            pad = 0.15
            shrink = 0.9
            rotation = 90

        cmap = mpl.cm.RdBu_r  # set the colormap to a divergent one

        indexes = np.arange(hessian.shape[0])

        x = [[i for i in indexes] for j in indexes]
        y = [[j for i in indexes] for j in indexes]

        lim = max(abs(min(hessian.flatten())), max(hessian.flatten()))

        im = ax.scatter(x, y, c=hessian.flatten(), marker='s',
                        cmap=cmap, vmin=-lim, vmax=lim)

        if cbar:
            cbar = plt.colorbar(im, ax=ax, format='%1.{}f'.format(deci),
                                orientation=orientation, pad=pad,
                                shrink=shrink)
            cbar.ax.tick_params(labelsize=ticks, rotation=rotation)
        return im

    def plot_ramachandran(self, step=1, marker_size_polar=5,
                          marker_size_rama=20, label_dots='Amino\nAcids'):
        """
        Shows the evolution of each phi-psi angle of each aminoacid in a polar
        and Ramachandran plot.

        Parameters
        ==========
        """
        rama_angles = self.rama_phi_psi([struct.atoms
                                         for struct in self.sith.structures])

        fig, axes = plt.subplots(3, 1)
        sp = StandardPlotter(fig=fig, ax=axes, figwidth=16,
                             figheight=22.5, ax_pref={'sci_not': False})
        scale = 0.05
        rs = np.arange(len(rama_angles))
        sp.set_polar(ax=0, r_ticks=rs[::step],
                     r_lims=[-rs[-1] * scale, rs[-1] * (1 + scale)])
        sp.set_polar(ax=1, r_ticks=rs[::step],
                     r_lims=[-rs[-1] * scale, rs[-1] * (1 + scale)])
        lines = []
        for j in range(len(rama_angles[0])):  # amino acids
            line, = sp.ax[0].plot(rama_angles[:, j][:, 0] * np.pi / 180, rs,
                                  '*',
                                  markersize=marker_size_polar,
                                  label=str(j + 1))
            lines.append(line)
            sp.ax[1].plot(rama_angles[:, j][:, 1] * np.pi / 180, rs, '*',
                          markersize=marker_size_polar)
            sp.ax[2].scatter(rama_angles[:, j][:, 0], rama_angles[:, j][:, 1],
                             s=marker_size_rama)

        vo = 5 / 8
        # adjust cartesian
        sp.add_space(borders=[[0, 0], [1, vo]],
                     axes=sp.ax[2])
        sp.spaces[-1].set_axis(borders=[[0, 0], [1, 1]],
                               spaces=(0.05, 0.1))
        sp.spaces[-1].set_axis(borders=[[0.16, 0.12], [0.98, 0.98]])
        # adjust polar
        sp.add_space(rows_cols=(1, 2),
                     borders=[[0, vo], [1, 1]],
                     axes=sp.ax[[0, 1]])
        sp.spaces[-1].set_axis(rows_cols=(1, 2),
                               borders=[[0.06, 0], [0.94, 0.94]],
                               spaces=(0.13, 0.1))

        leg = sp.spaces[-1].frame.legend(handles=lines,
                                         loc='upper center',
                                         bbox_to_anchor=[0.5, 0.4])
        leg.set_title(label_dots)

        sp.ax[0].set_title(r'$\phi$', fontsize=20)

        sp.ax[1].set_title(r'$\psi$', fontsize=20)

        sp.ax[0].set_rlabel_position(315)
        scale = 0.08
        sp.ax[0].set_rticks(rs[::step])
        sp.ax[0].set_ylim([-rs[-1] * scale, rs[-1] * (1 + scale)])

        sp.ax[1].set_rlabel_position(315)
        sp.ax[1].set_rticks(rs[::step])
        sp.ax[1].set_ylim([-rs[-1] * scale, rs[-1] * (1 + scale)])

        sp.ax[2].plot([0, 0], [-180, 180], color='gray')
        sp.ax[2].plot([-180, 180], [0, 0], color='gray')
        ticks = np.arange(-180, 180.1, 45, dtype=int)
        sp.ax[2].set_xticks(ticks)
        sp.ax[2].set_yticks(ticks)
        sp.ax[2].set_xlim([-180.1, 180.1])
        sp.ax[2].set_ylim([-180.1, 180.1])
        sp.ax[2].set_xlabel(r'$\phi$', fontsize=20)
        sp.ax[2].set_ylabel(r'$\psi$', fontsize=20)
        sp.ax[2].grid(True)
        sp.ax[2].tick_params(axis='both', labelsize=15)

        return sp
