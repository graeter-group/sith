from ase.calculators.gaussian import Gaussian
import numpy as np
from ase import Atom


class MoleculeSetter:
    """
    Tool to edit or extract information from a molecule.

    Parameters
    ==========
    atoms: ase.Atoms
        molecule to be analyzed.
    """
    def __init__(self, atoms):
        self.atoms = atoms

    def rot_x(self, angle):
        """
        Rotation matrix around x axis.

        Parameters
        ==========
        angle: float
            angle to rotate around the x axis in radians.

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[1, 0, 0],
                      [0, c, -s],
                      [0, s, c]])
        return R

    def rot_y(self, angle):
        """
        Retuns the rotation matrix around y axis.

        Parameters
        ==========
        angle: float
            angle to rotate around the y axis in radians.

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, 0, s],
                      [0, 1, 0],
                      [-s, 0, c]])
        return R

    def rot_z(self, angle):
        """
        Retuns the rotation matrix around z axis.

        Parameters
        ==========
        angle: float
            angle to rotate around the z axis in radians.

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector):
        """
        Apply the necessary rotations to set a vector aligned with positive x
        axis.

        Parameters
        ==========
        vector: array
            vector to be aligned.

        Return
        ======
        (numpy.array) [3float x 3float] Transformation matrix.
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
        Rotation around x axis to set a vector in the xy plane.

        Parameters
        ==========
        vector: array
            vector to be rotated to be in the xy plane.

        Return
        ======
        (numpy.array) [3float x 3float] Transformation matrix.
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1] / np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, trans, indexes=None, shift=None):
        """
        Apply a transformation to all vector positions of some atoms.

        Parameters
        ==========
        trans: array (3x3)
            transformation matrix to be applied to all atom positions.
        indexes: list or array. Defatult=None
            indexes of the atoms to apply the transformation. Default None
            that means the transformation is applied to the positions of all
            the atoms. 1-based numbering.
        shift: list. Default=None
            vector to apply a translation.

        Return
        ======
        (numpy.array)[#natoms x 3float] new xyz positions of the N atoms. It
            changes the positions of the internal atoms object.
        """
        if shift is None:
            shift = [0, 0, 0]
        if indexes is None:
            indexes = np.arange(len(self.atoms)) + 1

        new_positions = []
        for i, atom in enumerate(self.atoms):
            if (i + 1) in indexes:
                new_positions.append(np.dot(trans, atom.position) + shift)
            else:
                new_positions.append(atom.position + shift)
        self.atoms.set_positions(new_positions)

        return new_positions

    def xy_alignment(self, index1, index2, index3=None, center=None):
        """
        Transform the positions of the atoms such that the atoms of indexes1
        and index2 are aligned in the x axis. The atom with index3 would be in
        the xy plane in case to be given.

        Parameters
        ==========
        index1 and index2: int
            indexes of the atoms to be aligned with the x-axis. The positive
            direction of x would go from atom 1 to atom 2. 1-based numbering.
        index3: int. Default=None
            The atom with index 3 would be in the xy plane in case to be given.
            1-based numbering.
        Center: int. Default=None
            It must be index1 or index2, that means the atom with this index
            will be placed in the origin. In case center=None (default), the
            origin would be in the geometrical center between atoms with index1
            and index2.

        Return
        ======
        (numpy.array)[#natoms x 3float] new xyz positions of the N atoms. It
            changes the positions of the internal atoms object.
        """
        # Move the origin
        if center == index1 or center == index2:
            center = self.atoms[center - 1].position
        else:
            pos1 = self.atoms[index1 - 1].position
            pos2 = self.atoms[index2 - 1].position
            center = (pos1 + pos2) / 2

        self.atoms.set_positions(self.atoms.positions - center)
        # set index1 and index2 along x axis
        axis = (self.atoms[index2 - 1].position
                - self.atoms[index1 - 1].position)
        self.apply_trans(self.align_axis(axis))
        if index3 is not None:
            third = self.atoms[index3 - 1].position
            self.apply_trans(self.align_plane(third))
        return self.atoms

    def increase_distance(self, constraints, deltad):
        """
        Increases the distance between two atoms by moving them and keeping the
        rest of the atoms in the same place.

        Parameters
        ==========
        constraints:
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance. 1-based numbering.
        deltad: float
            amount to add to the distance between atoms.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        self.xy_alignment(constraints[0][0], constraints[0][1])

        d1norm = self.atoms.get_distance(constraints[0][0] - 1,
                                         constraints[0][1] - 1)
        self.atoms.set_distance(constraints[0][0] - 1, constraints[0][1] - 1,
                                d1norm + deltad)
        return self.atoms

    def increase_distance_with_constraints(self, constraints, deltad):
        """
        Takes a configuration and increases the distance between two atoms by
        moving those atoms and all those connected by constraints (all
        together) containing them and keeping the rest of the atoms in the same
        place.

        Parameters
        ==========
        constraints:
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance. 1-based numbering.
        deltad: float
            amount to add to the distance between atoms.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        self.xy_alignment(constraints[0][0], constraints[0][1])
        # Find all atoms connected by constraints to the first or second
        # atom of the first constraint:
        left = [constraints[0][0]]
        right = [constraints[0][1]]
        len_right = 0
        len_left = 0
        while ((len_left != len(left)) and (len_right != len(right))):
            len_left = len(left)
            len_right = len(right)
            for cons in constraints[1:]:
                if ((cons[0] in left) and (cons[1] not in left)):
                    left.append(cons[1])
                elif ((cons[1] in left) and (cons[0] not in left)):
                    left.append(cons[0])
                elif ((cons[0] in right) and (cons[1] not in right)):
                    right.append(cons[1])
                elif ((cons[1] in right) and (cons[0] not in right)):
                    right.append(cons[0])
        print("These will be the atoms moved to the left: ", left)
        print("These will be the atoms moved to the right: ", right)

        # Move to the left the atoms on the left and to the right the atoms
        # on the right, the rest of the atoms remain in the same place:
        new_positions = []
        for i, atom in enumerate(self.atoms):
            if i in np.array(right) - 1:
                deltapos = np.array([deltad / 2, 0, 0])
                new_positions.append(atom.position + deltapos)
            elif i in np.array(left) - 1:
                deltapos = np.array([-deltad / 2, 0, 0])
                new_positions.append(atom.position + deltapos)
            else:
                new_positions.append(atom.position)

        self.atoms.set_positions(new_positions)
        return self.atoms

    def scale_distance(self, constraints, deltad, index3=None):
        """
        Increase the distance between two atoms by aligning those atoms with
        the x axis and scaling the x-coordinate of all intermedia atoms.

        Parameters
        ==========
        constraints: list
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance. 1-based numbering.
        deltad: float
            amount to add to the distance between atoms.
        index3: int. Default=None
            The atom with index 3 would be in the xy plane in case to be given.
            1-based numbering.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        index1, index2 = constraints[0]
        d1norm = self.atoms.get_distance(index1 - 1, index2 - 1)
        # Move atom1 to the origin and rotate the molecule such that atom2 is
        # aligned with the +x axis:
        self.xy_alignment(index1, index2, center=index1)
        scale_factor = np.array([(d1norm + deltad) / d1norm, 1, 1])
        new_positions = [atom.position * scale_factor
                         for atom in self.atoms]
        self.atoms.set_positions(new_positions)
        self.xy_alignment(index1, index2, index3=index3)

        return new_positions

    def create_gaussian_input(self, out=None, charge=0, xc='bmk',
                              basis='6-31+g', mem='60GB'):
        """
        Creates a gaussian .com file without specifing the kind of calculus to
        run (optimization, ab-initio md, frequencies...). You would have to add
        it.

        Parameters
        ==========
        out: str. Defatult=None
            name of the gaussian file (.com) without extension. Default
            chemical formula.
        charge: int. Default=0
            charge of the molecule in electron units. Default 0.
        xc: str. Default='bmk'
            exchange correlation functional used in gaussian. Default bmk
        basis: str. Default='6-31+g'
            basis set used in gaussian. Default 6-31+g

        Return
        ======
        (ase.calculator.Gaussian) Calculator used to create the input.
        """
        if out is None:
            out = self.atoms.get_chemical_formula()

        calculator = Gaussian(mem=mem,
                              label=out,
                              chk=out,
                              xc=xc,
                              basis=basis,
                              mult=1,
                              charge=int(charge))

        calculator.write_input(self.atoms)

        return calculator


class Alignment:
    """
    Set of tools to apply the necessary transformations for aligning a the
    three first eigenvectors of PCA with the x, y, z axis, in that order.

    Parameters
    ==========
    atoms: ase.Atoms
        molecule to be aligned.
    """
    def __init__(self, atoms):
        self.atoms = atoms
        self.ms = MoleculeSetter(self.atoms)

    @staticmethod
    def pca_vectors(atoms, indexes='all'):
        """
        Find the eigen vectors of the Principal Component Analysis of the
        positis of the atoms.

        Parameters
        ==========
        atoms: ase.Atoms
            molecule to be analyzed.
        indexes: list. Default="all"
            subset of atoms to be included in the pca analysis. 1-based
            numbering.

        Return
        ======
        (array) principal components.
        """
        try:
            from sklearn.decomposition import PCA
        except ModuleNotFoundError:
            raise ModuleNotFoundError("This code requires sklearn")

        if indexes == 'all':
            indexes = np.arange(len(atoms)) + 1
        positions = atoms.positions[np.array(indexes) - 1]
        pca = PCA(n_components=3)
        pca.fit(positions)
        components = pca.components_
        return components

    @staticmethod
    def center_geo(atoms):
        """
        Moves the coordinates origin to the geometrical center atoms.

        Parameters
        ==========
        atoms: ase.Atoms
            Molecule to be moved.

        Return
        ======
        (ase.Atoms) Transformed atoms
        """
        ms = MoleculeSetter(atoms)
        geo_center = np.sum(atoms.positions, axis=0) / len(atoms)
        ms.apply_trans(np.identity(3), shift=-geo_center)

        return ms.atoms

    @staticmethod
    def align_with_components(atoms):
        """
        Aligns the coordinates x, y, z axis with the main vectors (in that
        order).

        Parameters
        ==========
        atoms: ase.Atoms
            Molecule to be aligned.

        Return
        ======
        (ase.Atoms) rotated and shifted molecule.
        """
        ms = MoleculeSetter(atoms)
        ms.atoms = Alignment.center_geo(ms.atoms)
        components = Alignment.pca_vectors(ms.atoms)

        # add dummy atoms
        ms.atoms += Atom('H', position=[0, 0, 0])
        for vec_comp in components:
            ms.atoms += Atom('H', position=vec_comp)

        # transform based on dummy atoms: center on the origin, first PCA
        # component with the positive x-axis, second PCA component in the
        # xy plane.
        ms.xy_alignment(-3, -2, -1)

        # remove dummy atoms
        ms.atoms = ms.atoms[:-4]

        return ms.atoms


class PCAMatcher:
    """
    Tools to align two molecules with their PCA main vectors.

    Parameters
    ==========
    reference: ase.Atoms
        reference of the molecule to be matched.
    to_compare: ase.Atoms
        molecule to be compared to the reference.
    """
    def __init__(self, reference, to_compare):
        assert len(reference) == len(to_compare), "The reference and the " + \
            "structure to compare have to have the same number of atoms"
        self.reference = reference.copy()
        self.to_compare = to_compare.copy()

    def align_molecules(self):
        """
        Align each molecule with their own main PCA vectors.

        Return
        ======
        (tuple) aligned reference and molecule to compare.
        """
        self.reference = Alignment.center_geo(self.reference)
        self.reference = Alignment.align_with_components(self.reference)
        self.to_compare = Alignment.center_geo(self.to_compare)
        self.to_compare = Alignment.align_with_components(self.to_compare)
        return self.reference, self.to_compare

    def matching(self):
        """
        Finds the order of the indices of the molecule which fit with the atoms
        of the reference. So, if the outcome of this function is "out", then
        PCAMatcher.to_compare[out[i]] == PCAMatcher.reference[i] should be
        true.

        Return
        ======
        (numpy.array) indices of the molecule that match the indices of the
        reference.

        Note
        ----
        align before using this function
        """
        self.align_molecules()
        n_atoms = len(self.reference)
        correspondence = -np.ones(n_atoms, dtype=int)

        for i in range(n_atoms):
            distances = np.linalg.norm(self.reference.positions[i]
                                       - self.to_compare.positions,
                                       axis=1)
            i_min = np.where(distances == min(distances))[0][0]

            if correspondence[i] == -1:
                correspondence[i] = int(i_min)

        if -1 in correspondence:
            print('Warning: check repetitions because the next atoms of the '
                  'reference did not get a correspondent in the test:',
                  np.where(correspondence == -1)[0])

        return correspondence
