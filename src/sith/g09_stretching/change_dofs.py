import numpy as np
from sith.utils.peptides import PepSetter
from sith.utils.miscellaneous import output_terminal
from ase.io import read, write
from pytest import approx
from sith.utils.tools import conf2pdb


def permute_atoms(atoms, indexes):
    """
    Permutes two atoms in an Atoms object.

    Parameters
    ==========
    atoms: ase.Atoms
        atoms object to be permuted.
    indexes: list
        list of two indexes to be permuted (1-based).

    Return
    ======
    (ase.Atoms) permuted atoms object.
    """

    indexes = list(indexes)
    indexes.sort()

    atoms = (atoms[:indexes[0] - 1]
            + atoms[indexes[1] - 1]
            + atoms[indexes[0]: indexes[1] - 1]
            + atoms[indexes[0] - 1]
            + atoms[indexes[1]:])

    return atoms


def extract_dofs(indexes, atoms):
    """
    Extracts the distance, angle and dihedral defined by the given indexes.

    Parameters
    ==========
    indexes: list
        list of four indexes (1-based) defining the distance (between the
        first two indices), angle (between the first three indices) and
        dihedral (the four indices) you want to compute.
    atoms: ase.Atoms
        atoms of a given configuration.

    Return
    ======
    (tuple) distance, angle and dihedral defined by the given indexes.
    """
    for i in range(len((indexes))):
        indexes[i] = indexes[i] - 1
    distance = atoms.get_distance(*indexes[:2])
    angle = atoms.get_angle(*indexes[:3])
    dihedral = atoms.get_dihedral(*indexes)
    # dihedral in range -180 to 180
    while dihedral > 180:
        dihedral -= 360

    return distance, angle, dihedral


def test_dofs(atoms, index, file):
    """"
    Tests that the distance, angle and dihedral defined by the given indexes
    correspond to those in the given file.

    Parameters
    ==========
    atoms: ase.Atoms
        atoms of a given configuration.
    indexes: list
        Indexes (1-based) of an atom in the z-matrix that contains the
        distance, angle and dihedral.
    file: str
        file where the distance, angle and dihedral are stored.

    Return
    ======
    (bool) True if the test passes. It raises an assertion error if the values
    do not
    """
    indexes = output_terminal(f"grep ',R{index},' {file} | " +
                            "awk -F ',' '{print $2, $4, $6}'").split()
    indexes = [int(i) for i in indexes]
    indexes = [index] + indexes

    dist, angl, dihe = extract_dofs(np.array(indexes), atoms)
    old_dist = float(output_terminal(f'sed -n "/R{indexes[0]}=/p" {file}'
                                    + ' | cut -d = -f 2', print_output=False))
    old_angl = float(output_terminal(f'sed -n "/A{indexes[0]}=/p" {file}'
                                    + ' | cut -d = -f 2', print_output=False))
    old_dihe = float(output_terminal(f'sed -n "/D{indexes[0]}=/p" {file}'
                                    + ' | cut -d = -f 2', print_output=False))

    assert old_dist == approx(dist, abs=3e-2), f"R{indexes[0]}({old_dist}) " +\
        f"does not correspond to the expected from the xyz file({dist})"
    assert old_angl == approx(angl, abs=5e-1), f"A{indexes[0]}({old_angl})" +\
        f"does not correspond to the expected from the xyz file({angl})"
    assert old_dihe == approx(dihe, abs=5e-1), f"D{indexes[0]}({old_dihe})" +\
        f"does not correspond to the expected from the xyz file({dihe})"

    return True


def def_line(indexes, element):
    """
    Creates a line defining an atom in a gaussian input comfile.

    Parameters
    ==========
    indexes: list
        list of four indexes (1-based) defining the distance (between the
        first two indices), angle (between the first three indices) and
        dihedral (the four indices).
    element: str
        element of the atom you want to define. E.g. 'C', 'H', 'O', etc.

    Return
    ======
    (str) line defining the atom in a gaussian input comfile.
    """
    a1, a2, a3, a4 = indexes
    return f'{element},{a2},R{a1},{a3},A{a1},{a4},D{a1},0'


def change_def(new_i, element, atoms, file):
    """
    Changes the definition of a given DOF and updates the comfile including the
    value of the DOFs.

    Parameters
    ==========
    new_i: list
        list of four indexes (1-based) defining the distance (between the
        first two indices), angle (between the first three indices) and
        dihedral (the four indices) you want to define.
    element: str
        element of the atom you want to define. E.g. 'C', 'H', 'O', etc.
    atoms: ase.Atoms
        atoms of a given configuration.
    file: str
        gaussian input comfile that you want to modify. Be sure to back up
        the original file.

    Return
    ======
    (None) However, it changes the comfile. Be sure to back up the original
    file.
    """
    # All the atoms used to define a new atom should already exist; they
    # should have lower index.
    for i in new_i[1:]:
        if new_i[0] < i:
            print("At least two atoms have to be swaped. You are trying to "
                + f"redefine {new_i}.")

    # Create new and old lines
    new_line = def_line(new_i, element)

    # Find variables
    tochange = new_i[0]
    new_i = np.array(new_i)
    dist, angl, dihe = extract_dofs(new_i, atoms)

    # Change value
    output_terminal(f'sed -i "/,R{tochange},/c\ {new_line}" {file}')
    output_terminal(f'sed -i "/R{tochange}=/c\ R{tochange}={dist}" {file}')
    output_terminal(f'sed -i "/A{tochange}=/c\ A{tochange}={angl}" {file}')
    output_terminal(f'sed -i "/D{tochange}=/c\ D{tochange}={dihe}" {file}')

# TODO: remove this function. it's deprecated
def find_CAC_NC(amino_info, atoms, i_pro):
    """
    Finds the C closest to the Ca atom and to the N atom. This is necessary
    because C could be in the same, the previous or in the next amino.

    Parameters
    ==========
    amino_info: dict
    dictionary with the amino info
    atoms: ase.Atoms
        atoms of a given configuration.
    i_pro: int
        index of the proline residue.

    Returns
    =======
    (tuple), index of C closest to Ca and index of the C closest to N. 1-based.
    """
    pre_amino = amino_info[i_pro - 1]
    amino = amino_info[i_pro]
    post_amino = amino_info[i_pro + 1]

    # obtain indices
    c_pre = pre_amino['C']
    c_cur = amino['C']
    c_pos = post_amino['C']
    ca_i = amino['CA']
    n_i = amino['N']

    # find Ca-C and N-C
    ca_dist = 100
    n_dist = 100
    cac = 0
    cn = 0
    for i in [c_pre, c_cur, c_pos]:
        # Ca-C
        d = atoms.get_distance(i - 1, ca_i - 1)
        if d < ca_dist:
            ca_dist = d.copy()
            cac = i
        # C-N
        d = atoms.get_distance(i - 1, n_i - 1)
        if d < n_dist:
            n_dist = d.copy()
            cn = i
    return cac, cn


def extract_proline_atoms(pep_set, i_pro):
    amino = pep_set.amino_info[i_pro]

    N_i = amino['N']
    Cd_i = amino['CD']
    Cg_i = amino['CG']
    Cb_i = amino['CB']
    Ca_i = amino['CA']
    HA_i = amino['HA']
    C_i = amino['C']
    O_i = amino['O']

    try:
        HB1_i = amino['1HB']
        HB2_i = amino['2HB']
        HG1_i = amino['1HG']
        HG2_i = amino['2HG']
        HD1_i = amino['1HD']
        HD2_i = amino['2HD']
    except KeyError:
        HB1_i = amino['HB1']
        HB2_i = amino['HB2']
        HG1_i = amino['HG1']
        HG2_i = amino['HG2']
        HD1_i = amino['HD1']
        HD2_i = amino['HD2']
    return N_i, Cd_i, Cg_i, Cb_i, Ca_i, HA_i, C_i, O_i, HB1_i, HB2_i, HG1_i, \
        HG2_i, HD1_i, HD2_i


def reorder_prolines_atoms(comfile, molecule, pdb_template, option):
    """
    Changes the definition of the degrees of freedom in prolines such that the
    missed bond is always different.

    Parameters
    ==========
    comfile: str
        gaussian input comfile that you want to modify.
    molecule: str
        xyz file representing the molecule.
    pdb_template: str
        pdb with the peptide information.
    option: str
        '1' misses Cgamma Cdelta.
        '2' misses Cbeta Cgamma.
        '3' misses Calpha Cbeta.
        '4' misses Cdelta N.

    Returns
    =======
    (None) However, it changes the comfile. Be sure to back up the original
    file.
    """
    pep_set = PepSetter(pdb_template)
    atoms = read(molecule)
    pros_i = np.where(np.array(list(pep_set.amino_name.values()))
                      == 'PRO')[0] + 1

    for i_pro in pros_i:
        pep_set = PepSetter(pdb_template)
        N_i, Cd_i, Cg_i, Cb_i, Ca_i, HA_i, C_i, O_i, HB1_i, HB2_i, HG1_i, \
            HG2_i, HD1_i, HD2_i = extract_proline_atoms(pep_set, i_pro)

        # test that xyz does does correspond to the com file
        test_dofs(atoms, Ca_i, comfile)
        test_dofs(atoms, Cb_i, comfile)
        test_dofs(atoms, Cg_i, comfile)
        test_dofs(atoms, Cd_i, comfile)
        test_dofs(atoms, N_i, comfile)
        test_dofs(atoms, HG1_i, comfile)
        test_dofs(atoms, HG2_i, comfile)

        if N_i < C_i:
            backbone = [N_i, Ca_i, HA_i, C_i, O_i]
        else:
            backbone = [C_i, O_i, N_i, Ca_i, HA_i]

        if option == '1' or option == '4':
            side = [Cb_i, HB1_i, HB2_i, Cg_i, HG1_i, HG2_i, Cd_i, HD1_i, HD2_i]
        elif option == '2':
            side = [Cb_i, HB1_i, HB2_i, Cd_i, HD1_i, HD2_i, Cg_i, HG1_i, HG2_i]
        elif option == '3':
            side = [Cd_i, HD1_i, HD2_i, Cg_i, HG1_i, HG2_i, Cb_i, HB1_i, HB2_i]
        
        reorder = np.array(backbone + side)

        new_atoms = conf2pdb(molecule, pdb_template, write_new_pdb=False)
        assert len(reorder) == max(reorder) - min(reorder) + 1, "The number of atoms in the " +\
            "pdb template and in the xyz file do not match."
        new_atoms = new_atoms[: min(reorder) - 1]  + new_atoms[reorder - 1] \
                    + new_atoms[max(reorder):]
        
        original_order = np.arange(min(reorder), max(reorder) + 1)
        print("finding permutations")
        permitations = swaps_to_transform(reorder,
                                          original_order)
        print("found permitations:", permitations)
        for swap in permitations:
            a = original_order[swap[0]]
            b = original_order[swap[1]]

            original_order[swap[0]], original_order[swap[1]] = b, a

            output_terminal(f"sith swap_atoms_in_com -a {a} " +\
                            f"-b {b} -f {comfile}")

        # save the new order
        pdb_template = molecule.replace('.xyz', '.pdb')
        write(pdb_template, new_atoms)

    return pdb_template

def swaps_to_transform(initial, final):
    """
    Return a list of swaps (i, j) that transform `initial` into `final`.
    Each swap is applied to `initial` in sequence.
    """
    arr = initial[:]          # make a copy
    pos = {v: i for i, v in enumerate(arr)}
    swaps = []

    for i in range(len(arr)):
        if arr[i] != final[i]:
            # Find where the desired element currently is
            j = pos[final[i]]

            # Record this swap
            swaps.append((i, j))

            # Perform the swap
            arr[i], arr[j] = arr[j], arr[i]

            # Update the position map
            pos[arr[j]] = j
            pos[arr[i]] = i

    return swaps

# add2executable
def change_prolines_dofs(comfile, molecule, pdb_template, option):
    """
    Changes the definition of the degrees of freedom in prolines such that the
    missed bond is always different.

    Parameters
    ==========
    comfile: str
        gaussian input comfile that you want to modify.
    molecule: str
        xyz file representing the molecule.
    pdb_template: str
        pdb with the peptide information.
    option: str
        '1' misses Cgamma Cdelta.
        '2' misses Cbeta Cgamma.
        '3' misses Calpha Cbeta.
        '4' misses Cdelta N.

    Returns
    =======
    (None) However, it changes the comfile. Be sure to back up the original
    file.

    Note
    ====
    This function only works when the prolines are not next to capping groups.
    """
    # first reorder the atoms in the prolines to avodid to define distances,
    # angles and dihedrals between atoms that are not yet defined in the
    # z-matrix
    pdb = reorder_prolines_atoms(comfile,
                                 molecule,
                                 pdb_template,
                                 option)
    # after reordering all the atoms in the prolines, change the dofs
    pep_set = PepSetter(pdb)
    atoms = read(pdb)
    write(molecule, atoms)
    pros_i = np.where(np.array(list(pep_set.amino_name.values()))
                      == 'PRO')[0] + 1

    for i_pro in pros_i:
        amino_prev = pep_set.amino_info[i_pro - 1]
        # extract the new definition of atoms
        N_i, Cd_i, Cg_i, Cb_i, Ca_i, HA_i, C_i, O_i, HB1_i, HB2_i, HG1_i, \
            HG2_i, HD1_i, HD2_i = extract_proline_atoms(pep_set, i_pro)
        
        if N_i < C_i:
            starting = 'N'
            previous = [amino_prev['C'], amino_prev['CA'], amino_prev['N']]
        else:
            starting = 'C'
            previous = [amino_prev['N'], amino_prev['CA'], amino_prev['C']]

        set_backbone(N_i, Ca_i, HA_i, C_i, O_i, previous, atoms, comfile,
                     starting)

        if option == '1':
            set1(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
                 HG2_i, HD1_i, HD2_i, atoms, comfile)
        if option == '2':
            set2(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
                 HG2_i, HD1_i, HD2_i, atoms, comfile)
        if option == '3':
            set3(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
                 HG2_i, HD1_i, HD2_i, atoms, comfile)
        if option == '4':
            set4(Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
                 HG2_i, HD1_i, HD2_i, atoms, comfile)
    
    return atoms


def set_backbone(N_i, Ca_i, HA_i, C_i, O_i, previous, atoms, comfile,
                 starting):
    if starting == 'N':
        change_def([N_i, previous[0], previous[1], previous[2]],
                   'N', atoms, comfile)
        change_def([Ca_i, N_i, previous[0], previous[1]],
                   'C', atoms, comfile)
        change_def([HA_i, Ca_i, N_i, previous[0]],
                   'H', atoms, comfile)
        change_def([C_i, Ca_i, N_i, previous[0]],
                   'C', atoms, comfile)
        change_def([O_i, C_i, Ca_i, N_i],
                   'O', atoms, comfile)
    elif starting == 'C':
        change_def([C_i, previous[0], previous[1], previous[2]],
                   'C', atoms, comfile)
        change_def([O_i, C_i, previous[0], previous[1]],
                   'O', atoms, comfile)
        change_def([Ca_i, C_i, previous[0], previous[1]],
                   'C', atoms, comfile)
        change_def([HA_i, Ca_i, C_i, previous[0]],
                   'H', atoms, comfile)
        change_def([N_i, Ca_i, C_i, previous[0]],
                   'N', atoms, comfile)

def set1(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
         HG2_i, HD1_i, HD2_i, atoms, comfile):
    change_def([Cb_i, Ca_i, C_i, O_i], 'C', atoms, comfile)
    change_def([HB1_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([HB2_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([Cg_i, Cb_i, Ca_i, C_i], 'C', atoms, comfile)
    change_def([HG1_i, Cg_i, Cb_i, Ca_i], 'H', atoms, comfile)
    change_def([HG2_i, Cg_i, Cb_i, Ca_i], 'H', atoms, comfile)
    change_def([Cd_i, N_i, Ca_i, C_i], 'C', atoms, comfile)
    change_def([HD1_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)
    change_def([HD2_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)


def set2(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
         HG2_i, HD1_i, HD2_i, atoms, comfile):
    change_def([Cb_i, Ca_i, C_i, O_i], 'C', atoms, comfile)
    change_def([HB1_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([HB2_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([Cd_i, N_i, Ca_i, C_i], 'C', atoms, comfile)
    change_def([HD1_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)
    change_def([HD2_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)
    change_def([Cg_i, Cd_i, N_i, Ca_i], 'C', atoms, comfile)
    change_def([HG1_i, Cg_i, Cd_i, N_i], 'H', atoms, comfile)
    change_def([HG2_i, Cg_i, Cd_i, N_i], 'H', atoms, comfile)

def set3(N_i, Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
         HG2_i, HD1_i, HD2_i, atoms, comfile):
    change_def([Cd_i, N_i, Ca_i, C_i], 'C', atoms, comfile)
    change_def([HD1_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)
    change_def([HD2_i, Cd_i, N_i, Ca_i], 'H', atoms, comfile)
    change_def([Cg_i, Cd_i, N_i, Ca_i], 'C', atoms, comfile)
    change_def([HG1_i, Cg_i, Cd_i, N_i], 'H', atoms, comfile)
    change_def([HG2_i, Cg_i, Cd_i, N_i], 'H', atoms, comfile)
    change_def([Cb_i, Cg_i, Cd_i, N_i], 'C', atoms, comfile)
    change_def([HB1_i, Cb_i, Cg_i, Cd_i], 'H', atoms, comfile)
    change_def([HB2_i, Cb_i, Cg_i, Cd_i], 'H', atoms, comfile)


def set4(Cd_i, Cg_i, Cb_i, Ca_i, C_i, O_i, HB1_i, HB2_i, HG1_i,
         HG2_i, HD1_i, HD2_i, atoms, comfile):
    change_def([Cb_i, Ca_i, C_i, O_i], 'C', atoms, comfile)
    change_def([HB1_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([HB2_i, Cb_i, Ca_i, C_i], 'H', atoms, comfile)
    change_def([Cg_i, Cb_i, Ca_i, C_i], 'C', atoms, comfile)
    change_def([HG1_i, Cg_i, Cb_i, Ca_i], 'H', atoms, comfile)
    change_def([HG2_i, Cg_i, Cb_i, Ca_i], 'H', atoms, comfile)
    change_def([Cd_i, Cg_i, Cb_i, Ca_i], 'C', atoms, comfile)
    change_def([HD1_i, Cd_i, Cg_i, Cb_i], 'H', atoms, comfile)
    change_def([HD2_i, Cd_i, Cg_i, Cb_i], 'H', atoms, comfile)