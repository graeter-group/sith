from typing import Union
import matplotlib as mpl
import numpy as np
from vmol.view import VMolecule
import vpython as vp
from sith.SITH import SITH
from sith.Utilities import color_distribution, create_colorbar


class EnergiesVMol(VMolecule):
    """
    Set of graphic tools to visualize the distribution of energies in the
    different degrees of freedom (lengths, angles, dihedrals). Wrap to VMol.

    Parameters
    ==========
    sith_info: SITH object.
        SITH object with the information of the system.
    dofs: list
        list of degrees of freedom defined according with 1-based indexing.
        It could also be 'all', 'bonds', 'angles' or 'dihedrals' to
        display all the corresponding DOFs.
    idef: int. default 0.
        index of the stretching to be displayed.
    alignment: list, tuple or np.ndarray. Default=None.
        indexes for the alignment of the molecule. See VMolecule.xy_alignment
        for more details.
    show_axis: bool. Default=False.
        True if you want to show the axis in the scene.
    background: list or vp.color. Default=vp.color.white.
        color of the background of the scene.
    \*\*kwargs: other arguments of VMolecule.
    """
    def __init__(self, sith_info: SITH, dofs: list,
                 idef: int = 0,
                 alignment: Union[list, tuple, np.ndarray] = None,
                 show_axis: bool = False,
                 background: Union[list, vp.color] = vp.color.white,
                 **kwargs):
        self._hook = False
        self.canvaskwargs = kwargs
        self.sith = sith_info
        self.energies = self.sith.dofs_energies
        dims = self.sith.dims
        self.nbonds = dims[1]
        self.nangles = dims[2]
        self.ndihedral = dims[3]

        atoms = [config.atoms for config in self.sith.structures]

        if idef < 0:
            assert abs(idef) <= len(atoms)
            self.idef = len(atoms) + idef
        else:
            self.idef = idef

        if 'height' in list(kwargs.keys()):
            height = kwargs['height']
        else:
            height = 500
        kwargs['height'] = height

        # matplotlib figure for colorbar
        self.fig = None
        self.kwargs_edofs = {'cmap': mpl.cm.get_cmap("Blues"),
                             'label': r"$\Delta$ Energy [Ha]",
                             'labelsize': 10,
                             'orientation': "vertical",
                             'div': 5,
                             'deci': 2,
                             'width': 700,
                             'height': 500}
        
        # The next fills the variable dofs with the definition of the degrees
        # of freedom
        # Create scene
        VMolecule.__init__(self, atoms,
                           show_axis=show_axis,
                           alignment=alignment,
                           frame=idef,
                           **kwargs)
        self._hook = False
        
        # create the list of dofs and color scale
        dofs = self._create_dofs_list(dofs)
        self.normalize, kwargs = self.create_figure(dofs, **kwargs)

        self.scene.background = self._asvector(background)
        self.traj_buttons()

        # show energies in dofs
        self.energies_some_dof(dofs, **self.kwargs_edofs)

        self._hook = True

    def _create_dofs_list(self, dofs: list) -> None:
        """
        Create the list of dofs to be displayed.

        Parameters
        ==========
        dofs: list of tuples or str.
            list of degrees of freedom defined according with 1-based indexing.
            It could also be 'all', 'bonds', 'angles' or 'dihedrals' to
            display all the corresponding DOFs.

        Return
        ======
        (list) list of lists with the dofs to be displayed according to SITH
        convention.
        """
        if 'all' in dofs:
            dofs = self.sith.dim_indices
        else:
            if 'bonds' in dofs:
                bonds = self.sith.dim_indices[:self.nbonds]
                dofs.extend(bonds)
                dofs.remove('bonds')
            elif 'angles' in dofs:
                angles = self.sith.dim_indices[self.nbonds:self.nbonds +
                                                   self.nangles]
                dofs.extend(angles)
                dofs.remove('angles')
            elif 'dihedrals' in dofs:
                dihedrals = self.sith.dim_indices[self.ndihedral:]
                dofs.extend(dihedrals)
                dofs.remove('dihedrals')

        return dofs

    
    def create_figure(self, dofs, **kwargs):
        """
        Creates the Color bar to be displayed at a side.

        Parameters
        ==========
        dofs: 
            DOFs to be displayed in the distribution.
        kwargs for change_def

        Return
        ======
        (list, kwargs) normalized energies, not used kwargs.
        """

        self.kwargs_edofs, kwargs = self.change_def(self.kwargs_edofs,
                                                    **kwargs)
        
        _, norm = color_distribution(self.sith,
                                     dofs,
                                     self.idef,
                                     self.kwargs_edofs['cmap'],
                                     absolute=True,
                                     div=self.kwargs_edofs['div'],
                                     decimals=self.kwargs_edofs['deci'])

        # Colorbar
        # Note that this colorbar 
        self.fig, _, _ = create_colorbar(norm, self.kwargs_edofs['label'],
                                        cmap=self.kwargs_edofs['cmap'],
                                        deci=self.kwargs_edofs['deci'],
                                        labelsize=self.kwargs_edofs['labelsize'],
                                        height=self.kwargs_edofs['height']/300,
                                        dpi=300)

        return norm, kwargs

    def energies_bonds(self, **kwargs) -> tuple:
        """
        Add the bonds with a color scale that represents the distribution of
        energy according to the JEDI method.

        Parameters
        ==========
        \*\*kwargs for EnergiesVMol.energies_some_dof

        Return
        ======
        (tuple) DOFs and their computed energies.
        """
        dofs = self.sith.dim_indices[:self.nbonds]
        out = self.energies_some_dof(dofs, **kwargs)

        return out

    def energies_angles(self, **kwargs) -> tuple:
        """
        Add the angles with a color scale that represents the distribution of
        energy according to the JEDI method.

        Parameters
        ==========
        \*\* kwargs for EnergiesVMol.energies_some_dof

        Return
        ======
        (tuple) DOFs and their computed energies.
        """
        dofs = self.sith.dim_indices[self.nbonds:self.nbonds +
                                                   self.nangles]
        out = self.energies_some_dof(dofs, **kwargs)
        return out

    def energies_dihedrals(self, **kwargs) -> tuple:
        """
        Add the dihedral angles with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        \*\* kwargs for EnergiesVMol.energies_some_dof

        Return
        ======
        (tuple) DOFs and their computed energies.
        """
        dofs = self.sith.dim_indices[self.nbonds + self.nangles:]
        out = self.energies_some_dof(dofs, **kwargs)
        return out

    def update_stretching(self, idef) -> None:
        """
        Function that is called when the frame is changed. Set the internal
        atribute self.kwargs_edofs to keep fixed the parameters to update the
        visualization of dofs.

        Parameters
        ==========
        idef: int
            index of the stretching to be updated to.

        Return
        ======
        (None)
        """
        self.idef = idef
        self.frame = idef
        self.atoms = self.trajectory[self.frame]

        for i, atom in enumerate(self.vatoms):
            atom.pos = vp.vector(*self.atoms[i].position)

        udofs = [dof.indices for dof in self.dofs.values()]
        if len(udofs) != 0:
            self.energies_some_dof(udofs, **self.kwargs_edofs)
        self.b_counter.text = f"   {self.idef}   "

    def energies_all_dof(self, **kwargs) -> tuple:
        """
        Add all DOF with a color scale that represents the distribution of
        energy according to the JEDI method.

        Parameters
        ==========
        \*\*kwargs for EnergiesVMol.energies_some_dof

        Return
        ======
        (tuple) DOFs and their computed energies.
        """
        dofs = self.sith.dim_indices
        return self.energies_some_dof(dofs, **kwargs)

    def change_def(self, def_dict: dict, **kwargs) -> tuple:
        """
        This functions change the values stored in a dictionary and removes
        each one of the arguments from the kwargs.

        Parameters
        ==========
        def_dict: dict
            dictionary with the default values.
        **kwargs: all the arguments you want to change.

        Return
        ======
        (dict, dict) modified dictionary withe the default values and set of
        kwargs without the used keys.
        """
        rem_keys = []
        for key, value in kwargs.items():
            if key in def_dict.keys():
                rem_keys.append(key)
                def_dict[key] = value

        for key in rem_keys:
            del kwargs[key]

        return def_dict, kwargs

    def energies_some_dof(self, dofs, **kwargs) -> tuple:
        """
        Add the bonds with a color scale that represents the distribution of
        energy according to the SITH method.

        Parameters
        ==========
        dofs: list of tuples.
            list of degrees of freedom defined according with 1-based indexing.
        normalize:
            normalization of colors acording to the range of energies and the
            colormap.
        \*\*kwargs of VMolecule.add_dof

        Return
        ======
        (list) VPython objects of the visual representation of DOFs.
        """
        self._hook = False
        self.kwargs_edofs, kwargs = self.change_def(self.kwargs_edofs,
                                                    **kwargs)
        cmap = self.kwargs_edofs['cmap']
        for dof in dofs:
            i = np.where(np.all(self.sith.dim_indices == dof, axis=1))[0]
            if len(i) != 1:
                raise ValueError(f"Dof {dof} not found or defined more than" +
                                 " once in the sith object.")
            color = cmap(self.normalize(self.energies[self.idef][i[0]]))[:3]
            self.add_dof(dof, color=color, **kwargs)
        inner_dofs = self.dofs
        self._hook = True

        return inner_dofs, kwargs

    def show_bonds_of_DOF(self, dof: list, unique: bool = False,
                          color: list = None) -> dict:
        """
        Show bonds of a specific dof.

        Parameters
        ==========
        dof: int.
            index in sith object that corresponds to the dof you want to show.
        unique: Bool. default False.. Default=False
            True if you want to remove all the other bonds and only keeping
            these ones.
        color: list[3(int)]. default R G B for angles, distances, dihedrals.
            color that you want to use in this dof.

        Return
        ======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        dof_indices = self.sith.dim_indices[dof]
        if color is None:
            colors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            color = colors[len(dof_indices) - 2]
        atoms1 = []
        atoms2 = []
        for i in range(len(dof_indices) - 1):
            atoms1.append(dof_indices[i])
            atoms2.append(dof_indices[i + 1])
        if unique:
            self.remove_all_bonds()

        return self.add_bonds(atoms1, atoms2, colors=color)

    def show_dof(self, dofs: list, **kwargs) -> dict:
        """
        Show specific degrees of freedom.

        Parameters
        ==========
        dofs: list of tuples.
            list of degrees of freedom defined according with 1-based indexing.
        kwargs for add_dof.

        Return
        ======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects

        Notes
        -----
        The color is not related with the JEDI method. It
        could be changed with the kwarg color=rgb list.
        """
        for dof in dofs:
            out = self.add_dof(dof, **kwargs)
        return out

    def show_bonds(self, **kwargs) -> None:
        """
        Show the bonds in the molecule.

        Parameters
        ==========
        kwargs for add_dof

        Return
        ======
        (None)

        Notes
        -----
        The color is not related with the JEDI method. It
        could be changed with the kwarg color=rgb list.
        """
        dofs = self.sith.dim_indices[:self.nbonds]
        self.show_dof(dofs, **kwargs)

    def traj_buttons(self) -> vp.button:
        """
        Create the buttons to move between stretching configurations.

        Return
        ======
        (vpython.button) Button that shows the index of the stretched
        configuration displayed.
        """
        def counter(b):
            pass

        def go_up(b, vmol=self):
            frame = vmol.idef
            if frame < len(vmol.trajectory) - 1:
                vmol.update_stretching(frame + 1)

        def go_down(b, vmol=self):
            frame = vmol.idef
            if frame > 0:
                vmol.update_stretching(frame - 1)

        # button to go down
        vp.button(text='\u25C4',
                  pos=self.scene.title_anchor,
                  bind=go_down)

        # counter box
        self.b_counter = vp.button(text=f"   {self.idef}   ",
                                   pos=self.scene.title_anchor,
                                   bind=counter)
        # button to go up
        vp.button(text='\u25BA',
                  pos=self.scene.title_anchor,
                  bind=go_up)

        return self.b_counter
