import matplotlib as mpl
from ipywidgets import HBox
import numpy as np
from sith.Utilities import color_distribution, create_colorbar


class MoleculeNGL:
    """
    Set of graphic tools to see the distribution of energies in the different
    degrees of freedom (lengths, angles, dihedrals). Wrap to NGLview.

    Parameters
    ==========
    atoms: ase.Atoms
        Atoms object to be visualized
    alignment: list[int]. Default=None
        list of three indexes corresponding to the
        indexes of the atoms in the xy plane. the first
        two atoms are set to the x axis.
    axis: bool. Default=False
        Show xyz axis.
    xsize: int. Default=500.
        Horizontal size of the visualization windows.
    ysize: int. Default=500.
        Vertical size of the visualization windows.
    """
    def __init__(self, atoms, alignment=None, axis=False,
                 xsize: int = 500, ysize: int = 500, n=5):
        try:
            from nglview import show_ase, show_asetraj
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("To use MoleculeNGL,"
                                      + " NGLview has to be installed. " + e)
        if type(atoms) is list:
            self.is_trajectory = True
            self.atoms = [config.copy() for config in atoms]
        else:
            self.is_trajectory = False
            self.atoms = atoms.copy()

        if alignment is not None:
            index1, index2, index3 = alignment
            if self.is_trajectory:
                [self.xy_alignment(config, index1, index2, index3)
                 for config in self.atoms]
            else:
                self.xy_alignment(self.atoms, index1, index2, index3)

        if self.is_trajectory:
            # Assume this is a trajectory or struct list
            self.view = show_asetraj(self.atoms, default=False)
            self.struct = self.atoms[0]
        else:
            # Assume this is just a single structure
            self.view = show_ase(self.atoms, default=False)
            self.struct = self.atoms.copy()

        self.view.add_spacefill()
        self.view.update_spacefill(radiusScale=0.15)
        self.view._remote_call('setSize', target='Widget',
                               args=['%dpx' % (xsize,), '%dpx' % (ysize,)])
        self.view.camera = 'orthographic'
        self.view.parameters = {"clipDist": 0}
        self.view.center()

        # The keys of the next dictionaries are the names of the DOFs and they
        # are defined such that the name is invariant to the order of the
        # indexes, such that the dof i-j((-k)-l) is the same as (l-(k-))j-i
        self.bonds = {}
        self.n = n
        self.angles = {}
        self.dihedrals = {}
        self.all_dofs_parameters = {}
        self.shape = self.view.shape
        self.box = self.view
        if axis:
            self.add_axis()

    def add_bond(self, atom1index, atom2index,
                 color=None, radius=0.1):
        """
        Add a bond between two atoms:
        atom1 and atom2

        Parameters
        ==========
        atom1index: int
            Indexes of the atoms to be connected using 1-based indexing.
        atom2index: int
            Indexes of the atoms to be connected using 1-based indexing.
        color: list. Default=gray([0.5, 0.5, 0.5])
            RGB triplet.
        radius: float. Default=0.1
            Radius of the bond.

        Return
        ======
        (Vpython object) Cylinder representing the just added bond.
        """
        if self.is_trajectory:
            atoms = self.atoms[self.view.frame]
        else:
            atoms = self.atoms
        if color is None:
            color = [0.5, 0.5, 0.5]

        indexes = [atom1index, atom2index]
        if atom1index > atom2index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)
        self.remove_bond(atom1index, atom2index)

        b = self.shape.add_cylinder(atoms[atom1index - 1].position,
                                    atoms[atom2index - 1].position,
                                    color,
                                    radius)

        self.bonds[name] = b
        self.all_dofs_parameters[name] = np.array([atom1index, atom2index,
                                                   0, 0])

        return self.bonds[name]

    def add_bonds(self, atoms1indexes, atoms2indexes, colors=None, radii=None):
        """
        Add a bond between each pair of atoms corresponding to
        two lists of atoms: atom1indexes and atom2indexes.

        Parameters
        ==========
        atoms1indexes: int
            Indexes of the atoms to be connected using 1-based indexing.
        atoms2indexes: int
            Indexes of the atoms to be connected using 1-based indexing.
        colors: list of color lists. Default all gray([0.5, 0.5, 0.5])
            RGB triplets for each of the bonds. It can be one a triplet
            in case of just one color in all bonds.
        radii: float or list of floats. Default 0.1
            radius of each bond.

        Return
        ======
        (list) Bonds in the system
        """

        if colors is None:
            colors = [0.5, 0.5, 0.5]

        if type(colors[0]) is not list:
            colors = [colors for i in range(len(atoms1indexes))]

        if radii is None:
            radii = 0.07

        if type(radii) is not list:
            radii = [radii for i in range(len(atoms1indexes))]

        assert len(atoms1indexes) == len(atoms2indexes), \
            "The number of atoms in both lists must be the same"
        assert len(atoms1indexes) == len(colors), \
            "The number of colors in must be the same as the number of atoms"
        assert len(atoms1indexes) == len(radii), \
            "The number of radii must be the same as the number of atoms"

        for i in range(len(atoms1indexes)):
            self.add_bond(atoms1indexes[i],
                          atoms2indexes[i],
                          colors[i],
                          radii[i])
        return self.bonds

    def remove_bond(self, atom1index, atom2index):
        """
        Remove a bond between two atoms:
        atoms1 and atoms2.

        Parameters
        ==========
        atom1index: int
            Indexes of the atoms that are connected. This bond
            will be removed.
        atom2index: int
            Indexes of the atoms that are connected. This bond
            will be removed.

        Return
        ======
        (list) Return the bonds in the system.
        """
        indexes = [atom1index, atom2index]
        indexes.sort()
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.bonds.keys():
            self.view.remove_component(self.bonds[name])
            del self.bonds[name]
            del self.all_dofs_parameters[name]
        return self.bonds

    def remove_bonds(self, atoms1indexes=None, atoms2indexes=None):
        """
        Remove several bonds in the plot between two list of atoms:
        atoms1 and atoms2.

        Parameters
        ==========
        atoms1indexes: list[int]. Default=None
            Indexes of the atoms that are connected.
        atoms2indexes: list[int]. Default=None
            Indexes of the atoms that are connected.

        Return
        ======
        (list) Return the bonds in the system.

        Note: if atoms2 is None, all bonds with atoms1 will me removed.
        If atoms1 and atoms2 are None, all bonds in the structure are
        removed.
        """

        if (atoms1indexes is None) and (atoms2indexes is None):
            for name in self.bonds.keys():
                self.view.remove_component(self.bonds[name])
            self.bonds.clear()
            return self.bonds

        elif (atoms1indexes is not None) and (atoms2indexes is None):
            to_remove = []
            for name in self.bonds.keys():
                for index in atoms1indexes:
                    if str(index) in name:
                        self.view.remove_component(self.bonds[name])
                        to_remove.append(name)
            for name in to_remove:
                del self.bonds[name]
                del self.all_dofs_parameters[name]
            return self.bonds

        else:
            assert len(atoms1indexes) == len(atoms2indexes), \
                "The number of atoms in both lists must be the same"
            [self.remove_bond(index1, index2)
             for index1, index2 in
             zip(atoms1indexes, atoms2indexes)]
            return self.bonds

    def remove_all_bonds(self):
        """
        Remove all bonds

        Return
        ======
        (list) Return the bonds in the system.
        """
        return self.remove_bonds()

    def add_arc(self, vertex, arcdots, color):
        """
        Add an arc using triangles.

        Parameters
        ==========
        vertex: array
            center of the arc
        arcdots: list of arrays
            vectors that define the points of the arc. These
            vectors must be defined respect the vertex.
        color: color format.
            color if the arc to be added

        Return
        ======
        (list) list of meshes forming the arc.
        """

        triangles = []
        for i in range(len(arcdots) - 1):
            vertexes = np.hstack((vertex,
                                  vertex + arcdots[i],
                                  vertex + arcdots[i + 1]))
            t = self.shape.add_mesh(vertexes, color)
            triangles.append(t)

        return triangles

    def add_angle(self, atom1index, atom2index, atom3index,
                  color=None, n=None):
        """
        Add an angle to between three atoms:
        atom1, atom2 and atom3
        - with the vertex in the atom2

        Parameters
        ==========
        atom1index: int
            Indexes of the first atom that defines the angle.
        atom2index: int
            Indexes of the second atom that defines the angle.
        atom3index: int
            Indexes of the third atom that defines the angle.
        color: color list. Default all gray([0.5, 0.5, 0.5])
            RGB triplet.
        n: int. Default=internal defined n
            number of intermedia points to add in the arc of
            the angle.

        Return
        ======
        (list) list of meshes forming the angle.
        """
        if n is None:
            n = self.n
        if color is None:
            color = [0.5, 0.5, 0.5]

        if self.is_trajectory:
            atoms = self.atoms[self.view.frame]
        else:
            atoms = self.atoms

        indexes = [atom1index, atom2index, atom3index]
        if atom1index > atom3index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)
        self.remove_angle(atom1index, atom2index, atom3index)
        self.angles[name] = []

        vertex = atoms[atom2index - 1].position
        side1 = atoms[atom1index - 1].position - vertex
        side2 = atoms[atom3index - 1].position - vertex
        lenside1 = np.linalg.norm(side1)
        lenside2 = np.linalg.norm(side2)
        lensides = min(lenside1, lenside2)
        side1 = 0.7 * lensides * side1 / lenside1
        side2 = 0.7 * lensides * side2 / lenside2

        arcdots = [side1, side2]
        color = color * 3

        new = self.intermedia_vectors(side1,
                                      side2,
                                      n)

        if n != 0:
            [arcdots.insert(1, vert) for vert in new[::-1]]

        self.angles[name] = self.add_arc(vertex, arcdots, color)
        self.all_dofs_parameters[name] = np.array([atom1index, atom2index,
                                                   atom3index, 0])

        return self.angles[name]

    def intermedia_vectors(self, a, b, n):
        """
        Define the intermedia arc dots between two vectors

        Parameters
        ==========
        a: array
            Vector of side a of the angle.
        b: array
            Vector of side b of the angle.
        n: int
             number of intermedia dots.

        Return
        ======
        (list) Return the intermedia vectors between two side vectors.
        """

        if n == 0:
            return []
        n += 1
        c = b - a
        lena = np.linalg.norm(a)
        lenb = np.linalg.norm(b)
        lenc = np.linalg.norm(c)
        lend = min(lena, lenb)

        theta_total = np.arccos(np.dot(a, b) / (lena * lenb))
        beta = np.arccos(np.dot(a, c) / (lena * lenc))
        intermedia = []
        for i in range(1, n):
            theta = i * theta_total / n
            gamma = beta - theta
            factor = (lena * np.sin(theta)) / (lenc * np.sin(gamma))
            dird = a + factor * c
            d = lend * dird / np.linalg.norm(dird)
            intermedia.append(d)
        return intermedia

    def remove_angle(self, atom1index, atom2index, atom3index):
        """
        Remove an angle if it exists

        Parameters
        ==========
        atom1index: int
            Indexes of the first atom that defines the angle.
        atom2index: int
            Indexes of the second atom that defines the angle.
        atom3index: int
            Indexes of the third atom that defines the angle.

        Return
        ======
        (list) Return all the angles in the display.
        """
        indexes = [atom1index, atom2index, atom3index]
        indexes.sort()
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.angles.keys():
            for triangle in self.angles[name]:
                self.view.remove_component(triangle)
            del self.angles[name]
            del self.all_dofs_parameters[name]

        return self.angles

    def remove_all_angles(self):
        """
        Remove all angles.

        Return
        ======
        (list) Return all the angles in the display.
        """
        names = self.angles.keys()

        for name in names:
            for triangle in self.angles[name]:
                self.view.remove_component(triangle)
        self.angles.clear()
        return self.angles

    def add_dihedral(self, atom1index, atom2index, atom3index,
                     atom4index, color=None, n=None):
        """
        Add an dihedral angle between four atoms:
        atom1, atom2, atom3 and atom4
        - with the vertex in the midle of the atom 2 and 3

        Parameters
        ==========
        atom1index: int
            Indexes of the first atom that defines the angle.
        atom2index: int
            Indexes of the second atom that defines the angle.
        atom3index: int
            Indexes of the third atom that defines the angle.
        atom4index: int
            Indexes of the fourth atom that defines the angle.
        color: color list. Default all gray([0.5, 0.5, 0.5])
            RGB triplet.
        n: int. Default 10
            number of intermedia points to add in the arc of
            the angle.

        Return
        ======
        (list) All the displayed dihedral angles.
        """
        if n is None:
            n = self.n

        if self.is_trajectory:
            atoms = self.atoms[self.view.frame]
        else:
            atoms = self.atoms
        if color is None:
            color = [0.5, 0.5, 0.5]

        indexes = [atom1index, atom2index, atom3index, atom4index]
        if atom1index > atom4index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        axis = (atoms[atom3index - 1].position
                - atoms[atom2index - 1].position)
        vertex = 0.5 * (atoms[atom3index - 1].position
                        + atoms[atom2index - 1].position)
        axis1 = (atoms[atom1index - 1].position
                 - atoms[atom2index - 1].position)
        axis2 = (atoms[atom4index - 1].position
                 - atoms[atom3index - 1].position)
        side1 = axis1 - axis * (np.dot(axis, axis1) / np.dot(axis, axis))
        side2 = axis2 - axis * (np.dot(axis, axis2) / np.dot(axis, axis))

        lenside1 = np.linalg.norm(side1)
        lenside2 = np.linalg.norm(side2)
        lensides = min(lenside1, lenside2)
        side1 = 0.7 * lensides * side1 / lenside1
        side2 = 0.7 * lensides * side2 / lenside2

        arcdots = [side1, side2]
        color = color * 3

        new = self.intermedia_vectors(side1,
                                      side2,
                                      n)

        if n != 0:
            [arcdots.insert(1, vert) for vert in new[::-1]]

        self.dihedrals[name] = self.add_arc(vertex, arcdots, color)
        self.all_dofs_parameters[name] = np.array([atom1index, atom2index,
                                                   atom3index, atom4index])

        return self.dihedrals[name]

    def remove_dihedral(self, atom1index, atom2index, atom3index, atom4index):
        """
        Remove the dihedral angle between 4 atoms:

        Parameters
        ==========
        atom1index: int
            Indexes of the first atom that defines the angle.
        atom2index: int
            Indexes of the second atom that defines the angle.
        atom3index: int
            Indexes of the third atom that defines the angle.
        atom4index: int
            Indexes of the fourth atom that defines the angle.

        Return
        ======
        (list) Return the dihedral angles.
        """
        indexes = [atom1index, atom2index, atom3index, atom4index]
        indexes.sort()
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.dihedrals.keys():
            for triangle in self.dihedrals[name]:
                self.view.remove_component(triangle)
            del self.dihedrals[name]
            del self.all_dofs_parameters[name]
        return self.dihedrals

    def remove_all_dihedrals(self):
        """
        Remove all dihedral angles.

        Return
        ======
        (list) all the dihedral angles
        """
        names = self.dihedrals.keys()

        for name in names:
            for triangle in self.dihedrals[name]:
                self.view.remove_component(triangle)
        self.dihedrals.clear()
        return self.dihedrals

    def add_dof(self, dof, color=None, **kwargs):
        """
        Add the degree of freedom to the molecule image

        Parameters
        ==========
        dof: tuple
            label of the degree of freedom using 1-based indexing.
        color: color format. Default=[0.5, 0.5, 0.5]
            Color of the DOF in the visual representation.
        kwargs:
            arguments for the function add_bond, add_angle or add_dihedral,
            according to the case.

        Return
        ======
        (list) list of DOFs of the same kind.

        Example
        =======
        dof=(1, 2) means a bond between atoms 1 and 2
        dof=(1, 2, 3) means an angle between atoms 1, 2 and 3
        dof=(1, 2, 3, 4) means a dihedral angle between atoms 1, 2, 3 and 4
        """

        if color is None:
            color = [0.5, 0.5, 0.5]

        types = ["bond", "angle", "dihedral"]
        type_dof = types[np.count_nonzero(dof) - 2]

        if type_dof == "bond":
            index1 = dof[0]
            index2 = dof[1]
            return self.add_bond(index1, index2, **kwargs)

        elif type_dof == "angle":
            index1 = dof[0]
            index2 = dof[1]
            index3 = dof[2]
            return self.add_angle(index1, index2, index3, **kwargs)

        elif type_dof == "dihedral":
            index1 = dof[0]
            index2 = dof[1]
            index3 = dof[2]
            index4 = dof[3]
            return self.add_dihedral(index1, index2, index3,
                                     index4, **kwargs)
        else:
            raise TypeError(f"{dof} is not an accepted degree of freedom.")

    def add_axis(self, length=1, radius=0.1):
        """
        Add xyz axis.

        Parameters
        ==========
        length: float. Default=1
            indicates the length of the axis in the visualization. Default=1
        radius: float. Default=0.1
            thickness of the xyz axis

        Return
        ======
        (dict) axis as dict with cylinder shapes of nglview.
        """
        self.axis = {}

        unit_vectors = np.array([[length, 0, 0],
                                 [0, length, 0],
                                 [0, 0, length]])
        for i in range(3):
            a = self.shape.add_cylinder([0, 0, 0],
                                        unit_vectors[i],
                                        unit_vectors[i] / length,
                                        radius)
            self.axis[str(i)] = a
        return self.axis

    def remove_axis(self):
        """
        remove xyz axis

        Return
        ======
        (dict) axis as dict with cylinder shapes of nglview.
        """
        for name in self.axis.keys():
            self.view.remove_component(self.axis[name])
        self.axis.clear()
        return self.axis

    def download_image(self, *args, **kwargs):
        """
        Download image.

        Parameters
        ==========
        args:
            args for the function nglview.view.download_image.
        kwargs:
            kwargs for the function nglview.view.download_image.

        Return
        ======
        (bool) True.
        """
        self.view.download_image(*args, **kwargs)
        return True

    def picked(self):
        """
        Call the function picked of nglview that checks the last clicked
        object and shows the information related to it.

        Return
        ======
        output of function nglview.view.picked
        """
        return self.view.picked

    # Alignment
    def rot_x(self, angle):
        """
        Returns the rotation matrix around x axis.

        Parameters
        ==========
        angle: float
            angle to rotate in radians.

        Return
        ======
        (np.array) 3x3 rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[1, 0, 0],
                      [0, c, -s],
                      [0, s, c]])
        return R

    def rot_y(self, angle):
        """
        Returns the rotation matrix around y axis.

        Parameters
        ==========
        angle: float
            angle to rotate in radians.

        Return
        ======
        (np.array) 3x3 rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, 0, s],
                      [0, 1, 0],
                      [-s, 0, c]])
        return R

    def rot_z(self, angle):
        """
        Returns the rotation matrix around z axis.

        Parameters
        ==========
        angle: float
            angle to rotate in radians.

        Return
        ======
        (np.array) 3x3 rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector):
        """
        Apply the necessary rotations to set a vector aligned with positive
        x axis.

        Parameters
        ==========
        vector: np.array
            vector to to be aligned.

        Return
        ======
        (np.array) 3x3 transformation matrix to align the vector with the x
        axis.
        """
        xyproj = vector.copy()
        xyproj[2] = 0
        phi = np.arcsin(vector[2] / np.linalg.norm(vector))
        theta = np.arccos(vector[0] / np.linalg.norm(xyproj))
        if vector[1] < 0:
            theta *= -1
        trans = np.dot(self.rot_y(phi), self.rot_z(-theta))
        return trans

    def align_plane(self, vector):
        """
        Rotation around x axis to set a vector in the xy plane

        Parameters
        ==========
        vector: np.array
            vector to be aligned.

        Return
        ======
        (np.array) 3x3 transformation matrix to align the vector with the x
        axis.
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1] / np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, atoms, trans, indexes=None):
        """
        Apply a transformation to all vector positions of the
        atoms object

        Parameters
        ==========
        atoms: ase.Atoms
            Atoms object to be transformed.
        trans:
            3x3 matrix containing the transformation.
        indexes: Default=None (namely, all atoms)
            indexes of atoms to apply the transformation.

        Return
        ======
        (np.array) array of new positions.
        """
        if indexes is None:
            indexes = list(range(len(atoms)))

        new_positions = []
        for i, atom in enumerate(atoms):
            if i in indexes:
                new_positions.append(np.dot(trans, atom.position))
            else:
                new_positions.append(atom.position)
        atoms.set_positions(new_positions)

        return new_positions

    def xy_alignment(self, atoms, index1, index2, index3):
        """
        Transforme the positions of the atoms such that
        the atoms of indexes 1 and 2 are aligned in the
        x axis

        Parameters
        ==========
        atoms: ase.Atoms
            Atoms to be transformed.
        index1: int
            first atom to go in the x axis.
        index2: int
            second atom to go to the x axis.
        index3: int
            atom to go in the xy plane.

        Return
        ======
        (np.array) numpy array with new positions.
        """
        # center
        center = (atoms[index1].position + atoms[index2].position) / 2
        atoms.set_positions(atoms.positions - center)
        axis = atoms[index2].position
        self.apply_trans(atoms, self.align_axis(axis))
        third = atoms[index3].position
        self.apply_trans(atoms, self.align_plane(third))

        return atoms.positions

    def show(self):
        """
        Show the molecule.

        Return
        ======
        (box)
        """
        return self.box


class EnergiesNGL(MoleculeNGL):
    """
    Set of tools to show a molecule and the distribution of energies in the
    different DOF.

    Parameters
    ==========
    sith_info:
        sith object or sith.utilities.ReadSummary object
    idef: int. Default='all'
        number of the deformation to be analized. Default=None, that means,
        all the structures are displayed as a trajectory.
    alignment: list. Default=None
        3 indexes to fix the correspondig atoms in the xy plane.
        The first atom is placed in the negative side of the x axis,
        the second atom is placed in the positive side of the x axis,
        and the third atom is placed in the positive side of the y axis.
    axis: bool
        add xyz axis
    background: color
        background color. Default: '#ffc'
    kwargs:
        kwargs of MoleculeNGL
    """
    def __init__(self, sith_info, idef='all', alignment=None, axis=False,
                 background='#FFFFFF', **kwargs):
        self.idef = idef
        self.sith = sith_info

        if idef == 'all':
            self.idef = None
            atoms = [config.atoms for config in self.sith.structures]
        else:
            self.idef = idef
            atoms = self.sith.structures[self.idef].atoms

        MoleculeNGL.__init__(self, atoms, alignment, axis, **kwargs)

        if self.idef is None:
            self.idef = self.view.frame

        self.view.background = background

        dims = self.sith.structures[0].dims
        self.nbonds = dims[1]
        self.nangles = dims[2]
        self.ndihedral = dims[3]
        self.kwargs_edofs = {'cmap': mpl.cm.get_cmap("Blues"),
                             'label': r"$\Delta$ Energy [Ha]",
                             'labelsize': 10,
                             'div': 5,
                             'deci': 2,
                             'width': 700,
                             'height': 500}
        self.view.observe(self.update_frame, names='frame')

    def energies_bonds(self, **kwargs):
        """
        Add the bonds with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        kwargs:
            of internal method energies_some_dof.

        Return
        ======
        (plt.figure) fig of energies_some_dof.
        """
        dofs = self.sith.structures[0].dim_indices[:self.nbonds]
        out = self.energies_some_dof(dofs, **kwargs)
        self.update_frame()
        return out

    def energies_angles(self, **kwargs):
        """
        Add the angles with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        kwargs:
            of internal method energies_some_dof.

        Return
        ======
        (plt.figure) fig of energies_some_dof.
        """
        dofs = self.sith.structures[0].dim_indices[self.nbonds:
                                                   self.nbonds + self.nangles]
        out = self.energies_some_dof(dofs, **kwargs)
        self.update_frame()
        return out

    def energies_dihedrals(self, **kwargs):
        """
        Add the dihedral angles with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        kwargs:
            of internal method energies_some_dof.

        Return
        ======
        (plt.figure) fig of energies_some_dof.
        """
        dofs = self.sith.structures[0].dim_indices[self.nbonds + self.nangles:]
        out = self.energies_some_dof(dofs, **kwargs)
        self.update_frame()
        return out

    def update_frame(self):
        """
        Function that is called when the frame is changed. Set the internal
        atribute self.kwargs_edofs to keep fixed the parameters to update the
        visualization of dofs.

        Return
        ======
        (bool) True
        """
        self.idef = self.view.frame

        dofs = list(self.all_dofs_parameters.values())
        if len(dofs) != 0:
            self.energies_some_dof(dofs, **self.kwargs_edofs)
        return True

    def energies_all_dof(self, **kwargs):
        """
        Add all DOF with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        kwargs:
            of internal method energies_some_dof.

        Return
        ======
        (plt.figure) fig of energies_some_dof.
        """
        dofs = self.sith.structures[0].dim_indices
        return self.energies_some_dof(dofs, **kwargs)

    def change_def(self, def_dict: dict, **kwargs) -> tuple:
        """
        This functions change the values stored in a dictionary and removes
        each one of the arguments from the kwargs.

        Parameters
        ==========
        def_dict: dict
            dictionary with the default values.
        kwargs:
            all the arguments you want to change.

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

    def energies_some_dof(self, dofs, **kwargs):
        """
        Add the bonds with a color scale that represents the
        distribution of energy according to the JEDI method.

        Parameters
        ==========
        dofs: list of tuples.
            list of degrees of freedom defined using 1-based indexing.
        kwargs:
            optional kwargs of internal method change_def

        Return
        ======
        (plt.figure) figure of the colorbar.
        """
        self.kwargs_edofs, kwargs = self.change_def(self.kwargs_edofs,
                                                    **kwargs)
        cmap = self.kwargs_edofs['cmap']
        label = self.kwargs_edofs['label']
        labelsize = self.kwargs_edofs['labelsize']
        div = self.kwargs_edofs['div']
        deci = self.kwargs_edofs['deci']
        width = self.kwargs_edofs['width']
        height = self.kwargs_edofs['height']

        energies, normalize = color_distribution(self.sith,
                                                 dofs,
                                                 self.idef,
                                                 cmap,
                                                 absolute=True,
                                                 div=div,
                                                 decimals=deci)

        for i, dof in enumerate(dofs):
            color = cmap(normalize(energies[self.idef][i]))[:3]
            self.add_dof(dof, color=color, **kwargs)

        # Insert colorbar in view
        self.view._remote_call("setSize",
                               targe="Widget",
                               args=[width, height])
        self.fig, self.ax, cbarwdg = create_colorbar(normalize,
                                                     label,
                                                     cmap,
                                                     deci,
                                                     labelsize,
                                                     height / 300)

        self.box = HBox([self.view, cbarwdg])

        return self.fig

    def show_bonds_of_DOF(self, dof, unique=False, color=None):
        """
        Show an specific dof.

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
        (list) Bonds in the system.
        """
        dof_indices = self.sith.structures[0].dim_indices[dof]
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

    def show_dof(self, dofs, **kwargs):
        """
        Show specific degrees of freedom.

        Parameters
        ==========
        dofs: list of tuples.
            list of degrees of freedom defined using 1-based indexing.
        kwargs for add_dof

        Return
        ======
        (None)

        Notes
        -----
        The color is not related with the JEDI method. It
        could be changed with the kwarg color=rgb list.
        """
        for dof in dofs:
            self.add_dof(dof, **kwargs)

    def show_bonds(self, **kwargs):
        """
        Show the bonds in the molecule.

        Parameters
        ==========
        kwargs for show_dof

        Return
        ======
        (None)

        Notes
        -----
        The color is not related with the JEDI method. It
        could be changed with the kwarg color=rgb list.
        """
        dofs = self.sith.structures[0].dim_indices[:self.nbonds]
        self.show_dof(dofs, **kwargs)
