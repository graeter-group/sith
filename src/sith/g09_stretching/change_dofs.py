import numpy as np
from sith.utils.peptides import PepSetter
from sith.utils.miscellaneous import output_terminal
from ase.io import read
from pytest import approx


def permute_atoms(atoms, indexes):
    indexes = list(indexes)
    indexes.sort()

    atoms = atoms[:indexes[0] - 1] + \
            atoms[indexes[1] - 1] + \
            atoms[indexes[0]: indexes[1] - 1] + \
            atoms[indexes[0] - 1] + \
            atoms[indexes[1]:]

    return atoms

def extract_dofs(indexes, atoms):
    for i in range(len((indexes))):
        indexes[i] = indexes[i] - 1
    distance = atoms.get_distance(*indexes[:2])
    angle = atoms.get_angle(*indexes[:3])
    dihedral = atoms.get_dihedral(*indexes)
    # dihedral in range -180 to 180
    while dihedral > 180:
        dihedral -= 360

    return distance, angle, dihedral


def test_dofs(atoms, indexes, file):
    dist, angl, dihe = extract_dofs(np.array(indexes) - 1, atoms)
    old_dist = float(output_terminal(f'sed -n "/R{indexes[0]}=/p" {file}' +
                                     ' | cut -d = -f 2', print_output=False))
    old_angl = float(output_terminal(f'sed -n "/A{indexes[0]}=/p" {file}' +
                                     ' | cut -d = -f 2', print_output=False))
    old_dihe = float(output_terminal(f'sed -n "/D{indexes[0]}=/p" {file}' +
                                     ' | cut -d = -f 2', print_output=False))

    assert old_dist == approx(dist, abs=3e-2), f"R{indexes[0]}({old_dist}) " +\
        f"does not correspond to the expected from the xyz file({dist})"
    assert old_angl == approx(angl, abs=5e-1), f"A{indexes[0]}({old_angl})" +\
        f"does not correspond to the expected from the xyz file({angl})"
    assert old_dihe == approx(dihe, abs=5e-1), f"D{indexes[0]}({old_dihe})" +\
        f"does not correspond to the expected from the xyz file({dihe})"

    return True


def def_line(indexes, element):
    a1, a2, a3, a4 = indexes
    print(f'new_line: {a1}, {a2}, {a3}, {a4}')
    return f'{element},{a2},R{a1},{a3},A{a1},{a4},D{a1},0'


def change_def(new_i, element, atoms, file):
    # All the atoms used to define a new atom should already exist; they
    # should have lower index.
    for i in new_i[1:]:
        if new_i[0] < i:
            print("At least two atoms have to be swaped. You are trying to " +
                  f"redefine {new_i}.")

    # Create new and old lines
    new_line = def_line(new_i, element)

    # Find variables
    tochange = new_i[0]
    new_i = np.array(new_i) - 1
    dist, angl, dihe = extract_dofs(new_i, atoms)

    # Change value
    output_terminal(f'sed -i "/,R{tochange},/c\ {new_line}" {file}')
    output_terminal(f'sed -i "/R{tochange}=/c\ R{tochange}={dist}" {file}')
    output_terminal(f'sed -i "/A{tochange}=/c\ A{tochange}={angl}" {file}')
    output_terminal(f'sed -i "/D{tochange}=/c\ D{tochange}={dihe}" {file}')


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
    (tuple), index of C closest to Ca and index of the C closest to N.
    """
    pre_amino = amino_info[i_pro - 1]
    amino = amino_info[i_pro]
    post_amino = amino_info[i_pro + 1]

    # obtain indices
    c_pre = pre_amino['C'] - 1
    c_cur = amino['C'] - 1
    c_pos = post_amino['C'] - 1
    ca_i = amino['CA'] - 1
    n_i = amino['N'] - 1

    # find Ca-C and N-C
    ca_dist = 100
    n_dist = 100
    cac = 0
    cn = 0
    for i in [c_pre, c_cur, c_pos]:
        # Ca-C
        d = atoms.get_distance(i, ca_i)
        if d < ca_dist:
            ca_dist = d.copy()
            cac = i + 1
        # C-N
        d = atoms.get_distance(i, n_i)
        if d < n_dist:
            n_dist = d.copy()
            cn = i + 1
    return cac, cn


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
    """
    pep_set = PepSetter(pdb_template)
    atoms = read(molecule)
    pros_i = np.where(np.array(list(pep_set.amino_name.values())) == 'PRO')[0] + 1

    for i_pro in pros_i:
        amino = pep_set.amino_info[i_pro]
        print(i_pro, amino)
        Ca_i = amino['CA']
        Cb_i = amino['CB']
        HB1_i = amino['2HB']
        HB2_i = amino['3HB']
        Cg_i = amino['CG']
        HG1_i = amino['2HG']
        HG2_i = amino['3HG']
        Cd_i = amino['CD']
        HD1_i = amino['2HD']
        HD2_i = amino['3HD']
        N_i = amino['N']
        CCA_i, CN_i = find_CAC_NC(pep_set.amino_info, atoms, i_pro)

        # test that xyz does does correspond to the com file
        test_dofs(atoms, np.array([Cg_i, Cb_i, Ca_i, N_i]), comfile)
        test_dofs(atoms, np.array([HG1_i, Cg_i, Cb_i, Ca_i]), comfile)
        test_dofs(atoms, np.array([HG2_i, Cg_i, Cb_i, HG1_i]), comfile)

        # option 1 (Default)
        set1(Ca_i, Cb_i, HB1_i, HB2_i, Cg_i, HG1_i, HG2_i, Cd_i, HD1_i,
             HD2_i, N_i, CCA_i, CN_i, atoms, comfile)

        if option == '2':
            # only changes Cgamma
            set2(Cg_i, Cd_i, N_i, Ca_i, HG1_i, HG2_i, atoms, comfile)
        
        if option == '3':
            # First change Cgamma
            atoms, Cg_i, Cd_i = set2(Cg_i, Cd_i, N_i, Ca_i, HG1_i, HG2_i, atoms, comfile)
            # and then change Cbeta
            set3(Cb_i, Cg_i, Cd_i, N_i, HB1_i, HB2_i, atoms, comfile)

        if option == '4':
            # Only changes Cdelta
            set4(Cd_i, Cg_i, Cb_i, Ca_i, HD1_i, HD2_i, atoms, comfile)


def set1(Ca_i, Cb_i, HB1_i, HB2_i, Cg_i, HG1_i, HG2_i, Cd_i, HD1_i,
            HD2_i, N_i, CCA_i, CN_i, atoms, comfile):
    print("Using set 1")

    # change betas dofs
    change_def(np.array([Cb_i, Ca_i, N_i, CN_i]),
               'C', atoms, comfile)
    change_def(np.array([HB1_i, Cb_i, Ca_i, N_i]),
               'H', atoms, comfile)
    change_def(np.array([HB2_i, Cb_i, Ca_i, N_i]),
               'H', atoms, comfile)

    # change gammas dofs
    change_def(np.array([Cg_i, Cb_i, Ca_i, N_i]),
               'C', atoms, comfile)
    change_def(np.array([HG1_i, Cg_i, Cb_i, Ca_i]),
               'H', atoms, comfile)
    change_def(np.array([HG2_i, Cg_i, Cb_i, Ca_i]),
               'H', atoms, comfile)

    # change deltas dofs
    change_def(np.array([Cd_i, N_i, Ca_i, CCA_i]),
               'C', atoms, comfile)
    change_def(np.array([HD1_i, Cd_i, N_i, Ca_i]),
               'H', atoms, comfile)
    change_def(np.array([HD2_i, Cd_i, N_i, Ca_i]),
               'H', atoms, comfile)


def set2(Cg_i, Cd_i, N_i, Ca_i, HG1_i, HG2_i, atoms, comfile):
    print("Using set 2")

    # change Cg dofs
    change_def(np.array([Cg_i, Cd_i, N_i, Ca_i]),
               'C', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HG1_i, Cg_i, Cd_i, N_i]),
               'H', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HG2_i, Cg_i, Cd_i, N_i]),
               'H', atoms, comfile)
 
    # swapping atoms Cg_i Cd_i
    print(f'swapping atoms {Cg_i} {Cd_i}')
    output_terminal(f'sith swap_atoms_in_com -a {Cg_i} -b {Cd_i} -f {comfile}')
    atoms = permute_atoms(atoms, [Cg_i, Cd_i])
    tmp  = Cg_i ; Cg_i = Cd_i ; Cd_i = tmp
    return atoms, Cg_i, Cd_i


def set3(Cb_i, Cg_i, Cd_i, N_i, HB1_i, HB2_i, atoms, comfile):
    print("Using set 3")

    # change Cg dofs
    change_def(np.array([Cb_i, Cg_i, Cd_i, N_i]),
               'C', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HB1_i, Cb_i, Cg_i, Cd_i]),
               'H', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HB2_i, Cb_i, Cg_i, Cd_i]),
               'H', atoms, comfile)
    
    print(f'swapping atoms {Cb_i} {Cg_i}')
    output_terminal(f'sith swap_atoms_in_com -a {Cb_i} -b {Cg_i} -f {comfile}')
    tmp  = Cb_i ; Cb_i = Cg_i ; Cg_i = tmp
    print(f'swapping atoms {Cg_i} {Cd_i}')
    output_terminal(f'sith swap_atoms_in_com -a {Cg_i} -b {Cd_i} -f {comfile}')


def set4(Cd_i, Cg_i, Cb_i, Ca_i, HD1_i, HD2_i, atoms, comfile):
    print("Using set 4")

    # change Cg dofs
    change_def(np.array([Cd_i, Cg_i, Cb_i, Ca_i]),
               'C', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HD1_i, Cd_i, Cg_i, Cb_i]),
               'H', atoms, comfile)

    # change Gg1 dofs
    change_def(np.array([HD2_i, Cd_i, Cg_i, Cb_i]),
               'H', atoms, comfile)