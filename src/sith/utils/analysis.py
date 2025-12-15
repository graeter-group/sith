# TODO: check indices in this file
import numpy as np
from sith import SITH
from sith.utils.peptides import PepSetter
from pathlib import Path
from typing import Union, Callable
from scipy.stats import f_oneway
import matplotlib.pyplot as plt
from PlottingTool import StandardPlotter


# region AnalyseData
def closest_indices(ref, value):
    """
    Find the indexes of the points after and before certain value

    Parameters
    ==========
    ref: numpy.array
        list contaning the data.
    value: float
        threshold value.

    Return
    ======
    (list) indexes of the next-lower and the next-upper elements in the list.
    """
    dif = ref - value
    multi = dif[1:] * dif[:-1]
    upp = np.where(multi <= 0)[0][-1] + 1
    low = upp - 1

    return [low, upp]


def cross_val(x, y, val):
    """
    Finds the crossing of a linearly interpolated curve with certain value in
    the y data (horizontal line).

    Parameters
    ==========
    x: numpy.array
        x-data of the curve.
    y: numpy.array
        y-data of the curve.
    val: float
        threshold value.

    Return
    ======
    (float) x value where the interpolated curve crosses the specified value.

    Note
    ====
    If you want to interpolate to a value of x (vertical line), just change x
    by y when you pass the arguments.
    """
    indices = closest_indices(y, val)
    x1, x2 = np.array(x)[indices]
    y1, y2 = np.array(y)[indices]

    x = (val - y1) * (x1 - x2) / (y1 - y2) + x1

    return x


def pep_per_x(xs, ys, value, dsa):
    """
    Separate the value of x in which y crosses a certain value, identifying it
    with an specific peptide.

    Parameters
    ==========
    xs: numpy.array
        x-data of the curve.
    ys: numpy.array
        y-data of the curve.
    value: float
        threshold value.
    dsa: DataSetAnalysis
        object containing the analysis.

    Return
    ======
    (dict) The name of the peptide as key and the x-value as dict value.
    """
    i = 0
    sep_info = {}
    for x, y in zip(xs, ys):
        try:
            x_cross = cross_val(x, y, value)
            pep = dsa.outcomes[i].name
            sep_info[pep] = x_cross
        except IndexError:
            print(f"{i} didn't even get the maximum energy, the maxener is ",
                  max(y))
        i += 1
    return sep_info


def reduce_data(xs, ys):
    """
    Reduces the data to the range of maximum minimum xs (highest of the lowest)
    and minumim maximum xs (lowest of the highest). Returns a new list of xs
    and ys where all the minum/maximum xs are the same and the ys of the
    minumum/maximum xs are found by interpolation. Also makes the
    transformation new_y -= new_y[:, 0]

    Parameters
    ==========
    xs: list of arrays
        set of xs data.
    ys: list of arrays
        set of xs data.

    Return
    ======
    (tuple) list of new xs and new ys

    Note
    ----
    if defmin/defmax is higher/lower than the real minimum/maximum, this code
    would not recognize the proper values. Be sure that defmin/defmax is
    lower/higher than the lowest/highest of xs.
    """
    minx = -np.inf
    maxx = np.inf

    for i, x in enumerate(xs):
        if x[0] > minx:
            minx = x[0]
        if x[-1] < maxx:
            maxx = x[-1]

    new_x = []
    new_y = []
    for i, x in enumerate(xs):
        miny = cross_val(ys[i], x, minx)
        maxy = cross_val(ys[i], x, maxx)
        condition = np.logical_and(x > minx, x < maxx)
        conca_x = np.append(x[condition], maxx)
        conca_y = np.append(ys[i][condition] - miny, maxy - miny)
        conca_x = np.insert(conca_x, 0, minx)
        conca_y = np.insert(conca_y, 0, 0)
        new_x.append(conca_x)
        new_y.append(conca_y)

    return new_x, new_y


def construct_neigh_dic(dsa, es):
    """
    Creates a dictionary with all the combination of two amino acids and stores
    in each combitation, a list of energies of the peptides with the first
    amino acid in the middle and the second in any side of the peptide.

    Parameters
    ==========
    dsa: DataSetAnalysis
        object that contains the information of the molecules in the dataset.
    es: list of arrays
        list of energies to be analyzed.

    Return
    ======
    (dict) keys: Amino combination, values: list of energies appearing in
    peptides with those combinations.

    Example
    =======
    # if energies per peptides are
    peptides = {'ABC': 1.5, 'GBA': 2.3}
    # neigh dictionary would be
    combi_neigh = {'BC': [1.5, 2.3]}
    """
    combi = {}
    for aa1 in np.unique(dsa.names.flatten()):
        for aa2 in np.unique(dsa.names.flatten()):
            combi[aa1 + aa2] = []

    for i, e in enumerate(es):
        pep = dsa.outcomes[i].name
        combi[pep[1] + pep[0]].append(e[-1])
        combi[pep[1] + pep[2]].append(e[-1])

    return combi


def separate_by_position(dsa, es, pos):
    """
    Stores the last values of "es" that correspond to an amino acid in a
    specific position.

    Parameters
    ==========
    dsa: DataSetAnalysis
        object that contains the information of the molecules in the dataset.
    es: list of arrays
        list of energies to be analyzed.
    pos: int
        position to be extracted.

    Return
    ======
    (dict) keys: amino acid in the asked position; values: list of energies to
    be analyzed.
    """
    assert len(es) == len(dsa.outcomes)

    energies = {}
    for aa in np.unique(dsa.names.flatten()):
        energies[aa] = []

    for i, e in enumerate(es):
        pep = dsa.outcomes[i].name
        energies[pep[pos]].append(e[-1])
    return energies


def mean_n_std(distri):
    """
    Extract the mean and standard deviation of a set of values associated to an
    amino acid in distri.

    Parameters
    ==========
    distri: dict
        keys amino acid, values list of values to obtain mean and standard
        deviation.

    Return
    ======
    (tuple) dictionary of mean, dictionary of standard deviation. Both
    dictionaries have the shape keys: amino acid; values: mean/std
    """
    mean_peps = {}
    std_peps = {}
    for aa in distri.keys():
        if len(distri[aa]) == 0:
            distri[aa] = None
            mean_peps[aa] = None
            std_peps[aa] = None
        else:
            mean_peps[aa] = np.mean(distri[aa])
            std_peps[aa] = np.std(distri[aa])
    return mean_peps, std_peps


def analysis_energy_dof(dsa, es):
    """
    Create a list of distributions per each one of the positions.

    Parameters
    ==========
    dsa: DataSetAnalysis
        object that contains the information of the molecules in the dataset.
    es: list of arrays
        list of energies to be analyzed.

    Return
    ======
    (duple) list of distribution of energies, mean and std per position.
    Namely, the first element of the first element of the output is the list of
    energies of the amino acid at the first position.
    """
    ener_per_pos1 = separate_by_position(dsa, es, 0)
    ener_per_pos2 = separate_by_position(dsa, es, 1)
    ener_per_pos3 = separate_by_position(dsa, es, 2)
    ener_per_pos = [ener_per_pos1,
                    ener_per_pos2,
                    ener_per_pos3]

    mean_per_pos = []
    stdr_per_pos = []
    for e_pos in ener_per_pos:
        mean, std = mean_n_std(e_pos)
        mean_per_pos.append(mean)
        stdr_per_pos.append(std)

    return ener_per_pos, mean_per_pos, stdr_per_pos


def statistics(ener_per_pos, pos, aminos):
    """
    Compute the p-value and F-statistic between distributions of the amino
    acids in the dataset for a given position.

    Parameters
    ==========
    ener_per_pos: list
        set of energy data per position, per amino acid.
    pos: int
        position to be analyzed.
    aminos: list
        list of amino acids in the data set.

    Return
    ======
    (tuple) arrays in shape of sqared matrices containing pval and fstat
    matrix.
    """
    pval_matrix = np.empty((len(aminos), len(aminos)))
    fstat_matrix = np.empty((len(aminos), len(aminos)))
    for i, aa1 in enumerate(aminos):
        for j, aa2 in enumerate(aminos):
            if (ener_per_pos[pos][aa2] is None) or \
               (ener_per_pos[pos][aa1] is None):
                pval_matrix[i][j] = None
                fstat_matrix[i][j] = None
                continue
            data1 = ener_per_pos[pos][aa1]
            data2 = ener_per_pos[pos][aa2]
            f_stat, p_value = f_oneway(data1, data2)
            pval_matrix[i][j] = p_value
            fstat_matrix[i][j] = f_stat

    return pval_matrix, fstat_matrix

# endregion


def dof_classificator_all(dofs_indexes, atoms_per_aminoacids):
    """
    Separates all degrees of freedom defined by atoms (all) of a residue.

    Parameters
    ==========
    dofs_indexes: list of duples
        sith.structures[n].dim_indices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (int), values: (array) [#DOFsPerResidue (int)]]
    Indexes of the degrees of freedom containing all atoms of each residue.

    Note
    ====
    atoms_per_aminoacids can be obtained from
    sith.utils.peptides.PepSetter.atom_indexes
    """
    list_aminos = {}
    for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
        list_aminos[i] = np.array([], dtype=int)
    for i in range(len(dofs_indexes)):
        for j in atoms_per_aminoacids.keys():
            if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).all():
                list_aminos[j] = np.append(list_aminos[j], i)
                break
    return list_aminos


def dof_classificator_one(dofs_indexes, atoms_per_aminoacids):
    """
    Separates all degrees of freedom defined by atoms (at least one) of a
    residue.

    Parameters
    ==========
    dofs_indexes: list of duples
        sith.structures[n].dim_indices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (str), values: (array) [#DOFsPerResidue (int)]]
    Indexes of the degrees of freedom containing at least one atom of each
    residue.

    Note
    ====
    atoms_per_aminoacids can be obtained from
    sith.utils.peptides.PepSetter.atom_indexes
    """
    list_aminos = {}
    for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
        list_aminos[i] = np.array([], dtype=int)
    for i in range(len(dofs_indexes)):
        for j in atoms_per_aminoacids.keys():
            if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).any():
                list_aminos[j] = np.append(list_aminos[j], i)
                break
    return list_aminos


def length_energy(sith, aminos_info, atom_types):
    """
    Return distances between two atom types in one amino acid and the
    energy associated with this DOF as the molecule is stretched.

    Parameters
    ==========
    sith: sith object
        sith object containing the distribution of energies. That implies to
        have the class variable 'energies' with the energies per deformed
        configuration and and per DOF.
    aminos_info: dic
        name of the atoms of one amino acid associated with the index.
    atom_types: str
        name of the atoms inside the aminoacid that will be studied,
        example ['CA', 'CB'].

    Return
    ======
    (list) [2 x #Def (float)] values of the DOF and energies associate with
    those DOFs per deformed configuration in the selected amino.

    Note
    ====
    aminos_info can be obtained from
    sith.utils.peptides.PepSetter.amino_info[n] where n is the selected amino
    acid.
    """
    defo = sith.structures[0]
    try:
        i_ric = defo.dim_indices.index((aminos_info[atom_types[0]],
                                       aminos_info[atom_types[1]]))
    except ValueError:
        i_ric = defo.dim_indices.index((aminos_info[atom_types[1]],
                                       aminos_info[atom_types[0]]))
    energies = sith.energies.T[i_ric]
    values_dof = []
    for defo in sith.structures:
        values_dof.append(defo.ric[i_ric])
    values_dof = np.array(values_dof)
    return [values_dof, energies]


# Deprecated
def le_same_aminoacids(sith, peptides_info, atom_types, kind_amino):
    """
    Return distances between two atom types in the same type of amino acid and
    the energy associated with these DOF as the molecule is stretched.

    Parameters
    ==========
    sith: sith object
        sith object containing the distribution of energies. That implies to
        have the class variable 'energies' with the energies per deformed
        configuration and and per DOF.
    peptides_info:
        object with the info of the peptide.
    atom_types: str
        name of the atoms inside the aminoacid that will be studied,
        example ['CA', 'CB'].
    kind_amino: (list) [(str)]
        name of the amino acids.

    Return
    ======
    (list) [#kindAmino x [2 x #Def (float)]] values of the DOF and energies
    associate with those DOFs per deformed configuration in the selected amino.

    Note
    ====
    peptides_info can be obtained from
    sith.utils.peptides.PepSetter
    """
    indexes = []
    for j, amino_name in enumerate(peptides_info.amino_name.values()):
        if amino_name in kind_amino:
            indexes.append(j + 1)
    all_le = []
    for index in indexes:
        values = length_energy(sith, peptides_info.amino_info[index],
                               atom_types)
        all_le.append(values)
    return all_le


class SithAnalysis:
    """
    Set of tools to anayze the data extracted from SITH.

    Parameters
    ==========
    sith: sith.SITH
       SITH object with all the information inside.
    pepinfo: dict
        Information of the peptide. Like amino acids, atoms, and so. It is an
        attribute of sith.utils.peptides.PepSetter.
    """
    def __init__(self, sith, pepinfo):
        self.sith = sith
        self.pep_info = pepinfo

    def le_dof_amino(self, a_names, aminos):
        """
        Value and the change of the energy of a DOF.

        Parameters
        ==========
        a_names: list
            atom names of the atoms forminf the DOF.
        aminos: list or int
            index of the residue from which each atoms belongs to. In case of
            only an integer, all the atoms are assumed to belong to the same
            residue.

        Return
        ======
        (tuple) array of dof values and change dof energies in the extension
        trajectory in respect to the optimized configuration.
        """
        # Note: aminos is actually residues
        if isinstance(aminos, int):
            # if all atoms belog to the same aminoacid
            aminos = [aminos for _ in a_names]
        else:
            assert len(aminos) == len(a_names)

        indexes = []
        for amino, atom in zip(aminos, a_names):
            indexes.append(self.pep_info.amino_info[amino][atom])

        dof = np.zeros(4, dtype=int)
        dof[-len(indexes):] = indexes
        dof_i = self.index_dof(dof)

        dof_value = self.sith.all_dofs[:, dof_i]
        dof_e = self.sith.dofs_energies[:, dof_i]
        dof_e -= dof_e[0]

        return dof_value, dof_e

    def index_dof(self, target: tuple):
        """
        Search the index of a specific dof.

        Parameters
        ==========
        target: np.ndarray
            degree of freedom.

        Return
        ======
        (int) index of the DOF.
        """
        for i, dof in enumerate(self.sith.dim_indices):
            dof_wo_0 = dof[np.nonzero(dof)[0]]
            target_wo_0 = target[np.nonzero(target)[0]]

            if len(dof_wo_0) != len(target_wo_0):
                continue

            if (dof_wo_0 == target_wo_0).all() or \
               (dof_wo_0 == target_wo_0[::-1]).all():
                return i
        raise ValueError(f"Non-found dof ({target}).")


def set_hes_from_ref(geo_ref, sith_tar, structure):
    """
    Set the hessian in in a target sith taken from a geometry of reference.

    Parameters
    ==========
    geo_ref: sith.Utilities.Geometry
        sith object that contains the atribute you want to redefine.
    sith_tar: sith.SITH
        sith object that will change its property.
    structure: int
        index of the deformed structure to set the hessian.

    Return
    ======
    (sith.SITH) returns the sith_tar with the hessian in the defined structure.

    Note: All the sith.sith.structures are Geometry objects with all the
    information of the structure.
    """
    for dof in sith_tar.dim_indices:
        test = dof[dof != 0]
        check2 = np.concatenate((test[::-1], np.zeros(4 - len(test),
                                                      dtype=int)))
        if not (np.all(geo_ref.dim_indices == dof, axis=1).any()
                or np.all(geo_ref.dim_indices == check2, axis=1).any()):
            raise ('this dof does not exist: ', dof)

    order = []
    for dof in sith_tar.dim_indices:
        test = dof[dof != 0]
        check2 = np.concatenate((test[::-1], np.zeros(4 - len(test),
                                                      dtype=int)))
        try:
            index = np.where(np.all(geo_ref.dim_indices == dof, axis=1))[0][0]
        except IndexError:
            index = np.where(np.all(geo_ref.dim_indices == check2,
                                    axis=1))[0][0]
        order.append(index)

    geo_ref.hessian = geo_ref.hessian[order]
    geo_ref.hessian = geo_ref.hessian[:, order]

    sith_tar.structures[structure].hessian = geo_ref.hessian

    return sith_tar


class DataSetAnalysis:
    """
    Creates the objects of to do the analysis of a set of peptides in a
    data set.

    Parameters
    ==========
    inner_steps: function. Default=None
        function with the operations to apply to every sith object -its
        unique argument- after initialized. Note that it has to have the
        analysis you want to apply. For example,
    data_dir: str. Default='./'
        directory containing the directories of the molecules.
    subdir: str. Default=''
        subdirectory inside each molecule directory containing the files to be
        analyzed by sith.
    exclude_prolines: bool. Default=True
        if True, neglects all the peptides with a proline in their sequence. In
        other words, it ignores all the directories with a 'P' in their name.
    exclude: list. Default=None
        list of molecules to exclude from the analysis.
    pdb_pattern: str. Default='stretched00'
        pattern to identify the pdb file in each molecule directory.
    struc_pattern: Default='forces<all>.fchk'
        pattern to identify the files to be analyzed by sith in each molecule
        directory.
    kwargs:
        for initializing sith


    .. note::

        The function inner_steps has to have as unique argument a sith object
        and return the same sith object after applying the analysis you want.
        For example, if you want to remove the first and last structure and
        then apply the sith_analysis, you can define the function as:

        .. code-block:: python

            def inner_steps(sith):
                sith.rem_first_last(from_last_minimum=True)
                sith.sith_analysis()
                return sith
    """

    def __init__(self, inner_steps: Callable, data_dir: str = './',
                 subdir='',
                 exclude_prolines=True,
                 exclude=None,
                 pdb_pattern='stretched00',
                 struc_pattern='force*.fchk',
                 alternative_subdir=None,
                 **kwargs):
        if isinstance(data_dir, (Path, str)):
            path = Path(data_dir)
        assert path.is_dir(), f"{data_dir} does not exist."

        peptides = list(path.glob('*/'))
        peptides.sort()
        assert len(peptides) != 0, f"There are not directories in {data_dir}"

        self.pep_infos = []
        self.names = []
        self.outcomes = []
        self.analysis = []

        print("Log SITH analysis:\n")

        prolines = []
        errors = []
        eff_exclu = []
        if exclude is None:
            exclude = []

        n_pep = 0
        for pep in peptides:
            if ((exclude_prolines) and ('P' in pep.stem)):
                prolines.append(pep.stem)
                continue
            if (pep.stem in exclude):
                eff_exclu.append(pep.stem)
                continue
            if pep.is_file():
                continue

            assert pep.is_dir(), f"{pep} does not exist."

            pdb = list(pep.glob(f'*{pdb_pattern}*.pdb'))
            if len(pdb) == 0:
                raise FileNotFoundError(f"Not '*{pdb_pattern}*.pdb' found in"
                                        f"{str(pep)}")
            elif len(pdb) > 1:
                raise ValueError("There are more than one pdb files with the "
                                 f"pattern {str(pep)}/*{pdb_pattern}*.pdb")

            self.pep_infos.append(PepSetter(pdb[0]))

            struc = pep / subdir
            structure_files = list(struc.glob(f'*{struc_pattern}*'))
            structure_files.sort()

            try:
                sith = SITH(inputfiles=structure_files, **kwargs)
                sith = inner_steps(sith)
                sith.name = pep.stem
                self.outcomes.append(sith)
                self.analysis.append(SithAnalysis(self.outcomes[-1],
                                                  self.pep_infos[-1]))
                if not self.test_continuity(self.outcomes[-1],
                                            self.analysis[-1]):
                    self.outcomes.pop(-1)
                    self.analysis.pop(-1)
                    raise ValueError("Discontinuity detected")

            except Exception as e:
                print(f"\nError in {str(pep)}:", e)
                try:
                    if alternative_subdir is None:
                        raise ValueError("No alternative subdir defined.")
                    struc = pep / alternative_subdir
                    structure_files = list(struc.glob(f'*{struc_pattern}*'))
                    structure_files.sort()
                    sith = SITH(inputfiles=structure_files, **kwargs)
                    sith = inner_steps(sith)
                    sith.name = pep.stem
                    self.outcomes.append(sith)
                    self.analysis.append(SithAnalysis(self.outcomes[-1],
                                                      self.pep_infos[-1]))
                    if not self.test_continuity(self.outcomes[-1],
                                                self.analysis[-1]):
                        self.outcomes.pop(-1)
                        self.analysis.pop(-1)
                        raise ValueError("Discontinuity detected")
                except Exception as e2:
                    print(f"\nError also in subdir:", e2)
                    self.pep_infos.pop(-1)
                    errors.append(pep.stem)
            assert len(self.pep_infos[-1].atoms)   == len(self.outcomes[-1].structures[0].atoms), f"Atom count mismatch for {pep.stem}"

            if n_pep % 20 == 0:
                print()
            n_pep += 1

            print(pep.stem + ' ', end='')

        print(f"\n--- A total of {len(self.outcomes)} peptides where added "
              "to the analysis")

        if len(prolines) > 0:
            print("--- The next peptides were neglected because they have a "
                  "proline at least:")
            [print(pep + ' ', end='') for pep in prolines]

        if len(errors) > 0:
            print("\n--- The next peptides did not worked for some reason. "
                  "Check them individually:")
            [print(pep + ' ', end='') for pep in errors]
        if len(eff_exclu) > 0:
            print("--- The next peptides were neglected because they have a "
                  "proline at least:")
            [print(pep + ' ', end='') for pep in eff_exclu]

        self.test()

    def test_continuity(self, sith, analysis):
        # atoms
        ca = analysis.pep_info.amino_info[3]['CA']
        c = analysis.pep_info.amino_info[3]['C']
        n = analysis.pep_info.amino_info[3]['N']

        # angle and distances
        ang = analysis.index_dof(np.array([c, ca, n]))
        cac = analysis.index_dof(np.array([c, ca]))
        can = analysis.index_dof(np.array([ca, n]))

        delta_angle = max(abs(sith.delta_q[:, ang].flatten() * 180 / np.pi))
        delta_dista = max(abs(sith.delta_q[:, [cac, can]].flatten()))
        if delta_angle > 10 or delta_dista > 0.2:
            return False
        return True

    def test(self):
        """
        check that the first sith has the basic variables. It assumes that
        the rest also has it

        Return
        ======
        (None)
        """
        assert isinstance(self.outcomes[0].structures_scf_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dims, np.ndarray)
        assert isinstance(self.outcomes[0].structure_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dofs_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dim_indices, np.ndarray)
        assert isinstance(self.outcomes[0].structures, list)

    def _populate_variables(self):
        """
        This method creates the attributes used fot the analysis extracting the
        values from the sith_objects.

        Return
        ======
        (DataSetAnalysis) This same object.
        """
        self.all_dft_energies = []
        self.all_edm_energies = []  # energy distribution method
        for ed_obj in self.outcomes:
            self.all_dft_energies.append(ed_obj.structures_scf_energies)
            self.all_edm_energies.append(ed_obj.structure_energies)

        self.all_dft_energies = np.array(self.all_dft_energies)
        self.all_edm_energies = np.array(self.all_edm_energies)

        return self

    def dof_vs_energy(self, aa_names, naas=3):
        """
        Extract the values of the dofs and the associated dofs energies
        for a set all the molecules of data set. Check
        sith.sith.analysis.SithAnalysis.le_dof_amino

        Parameters
        ==========
        aa_names: list
            Atom names forming the DOF.
        naas: list or int. Default=3
            index of the residue from which each atoms belongs to. In case of
            only an integer, all the atoms are assumed to belong to the same
            residue.

        Return
        ======
        (tuple) list of arrays  with dof values and dof energies in the
        extension trajectory of each molecule in the dataset.

        Note
        ----
        This function returns the absolute value of the DOF, to obtain the
        same but computing the change in the DOF in respect to the optimized
        configuration of see
        sith.sith.analysis.DataSetAnalysis.dofs_vs_energy_rescaled
        """
        xs = []
        es = []
        for an in self.analysis:
            x, e = an.le_dof_amino(aa_names, naas)
            xs.append(x)
            es.append(e)

        return xs, es

    def dof_vs_energy_rescaled(self, aa_names, naas=3):
        """
        Extract the values of the delta dofs in respect to the optimized
        configuration and the associated dofs energies for a set all the
        molecules of data set. Check
        sith.sith.analysis.SithAnalysis.le_dof_amino

        Parameters
        ==========
        aa_names: list
            Atom names forming the DOF.
        naas: list or int. Default=3
            index of the residue from which each atoms belongs to. In case of
            only an integer, all the atoms are assumed to belong to the same
            residue.

        Return
        ======
        (tuple) list of arrays with delta dof values and dof energies in the
        extension trajectory in respect to the optimized configuration of each
        molecule in the dataset.
        """
        xs = []
        es = []
        for an in self.analysis:
            x, e = an.le_dof_amino(aa_names, naas)
            xs.append(x - x[0])
            es.append(e)

        return xs, es

    def plot_DFT_ener(self, ax: plt.Axes = None, sp=None,
                      lw=1, ms=1, **kwargs):
        """
        Plots the DFT energy of the extension process

        Parameters
        ==========
        ax: plt.Axes Default=None
            Axes to plot in. If Nonne, a new axis is created.
        sp: StandardPlotter. Default=None
            sith.plotters.StandardPlotter
        lw: float. Default=1
            line width.
        ms: float. Default=1
            marker size.
        kwargs for StandardPlotters.

        Return
        ======
        (StandardPlotter, list, list) Standard plotter, all dofs values,
        all energy values used for the plot.
        """
        if 'ax_pref' in kwargs.keys():
            setter = kwargs['ax_pref']
            del kwargs['ax_pref']
        else:
            setter = {}

        if sp is None:
            sp = StandardPlotter(**kwargs)
        if ax is None:
            ax = sp.ax[0]
        sp.axis_setter(ax=ax,
                       xlabel=f'Stretched Structure',
                       ylabel='DFT Energy[Ha]',
                       **setter)

        xs = []
        ys = []
        for i, sith in enumerate(self.outcomes):
            y = sith.structures_scf_energies
            index1 = self.pep_infos[i].amino_info[1]['CH3']
            index2 = self.pep_infos[i].amino_info[5]['CH3']
            x = []
            for struc in sith.structures:
                try:
                    dist = struc.atoms.get_distance(index1 - 1, index2 - 1)
                except IndexError:
                    print(f"Error getting distance for {sith.name} setting distance to 0.0")
                    continue
                x.append(dist)
            sp.plot_data(x, y, ax=ax, lw=lw, markersize=ms)
            xs.append(x)
            ys.append(y)

        return sp, xs, ys

    def amino_freq(self):
        """
        Extract the name of the peptides in one-letter-convention for
        the amino acids.

        Return
        ======
        (numpy.array) set of names.
        """
        self.names = []
        for sith in self.outcomes:
            self.names.append(list(sith.name))
        self.names = np.array(self.names)

        return self.names
