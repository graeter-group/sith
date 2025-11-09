from sith.utils.molecules import Alignment, MoleculeSetter
from ase.geometry.analysis import Analysis
from ase.io import read, write
import numpy as np
import glob
from ase.constraints import FixAtoms


# add2executable
def shake_except(xyz_file, file_cons, modify_input=True, stdev=0.05):
    """
    Modify a xyz file by adding random noise to the position of the atoms
    except for those listed in the first line of a file.

    Parameters
    ==========
    xyz_file: int
        original xyz file to be modified.
    file_cons:
        file containing the atoms (1-based indices) to avoid to be modified by
        the shaking.
    modify_input: bool. Default=True
        True to modify the input file.
    scale: float. Default=0.05
        max magnitud of the noise to be added in each coordinate of each atom
        in Angstrom.
    stdev: float. Default=0.05
        standard deviation of the normal distribution used to generate the
        random noise.

    Return
    ======
    (ase.Atoms) new object atoms with the modified positions.
    """
    atoms = read(xyz_file)
    cons = np.loadtxt(file_cons, usecols=[0, 1], dtype=int) - 1
    if type(cons[0]) is np.int64:
        cons = np.array([list(cons)])

    c = FixAtoms(indices=cons[0])
    atoms.set_constraint(c)
    atoms.rattle(stdev=stdev, rng=np.random)

    manipulator = MoleculeSetter(atoms)
    manipulator.xy_alignment(cons[0][0], cons[0][1])

    if modify_input:
        write(xyz_file, manipulator.atoms)

    return atoms


# add2executable
def change_distance(inp, out, file_cons, deltad, charge, method, **kwargs):
    """
    Takes a configuration and increase the distance between two atoms. With the
    new structure, it creates a gaussian file without specifing the kind of
    calculation to run (optimization, ab-initio md, frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    file_cons: str
        file with the constraints sorted in the first two columns. The two
        atoms in the first row are the selected to change the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in electron charge units to create the gaussian input file. The
        multiplicity is assumed to be one.
    method: str
        method defined in sith.ase_utils.tools. So far, the methods already
        implemented are: 'scale_distance', 'increase_distance',
        'increase_distance_with_constraints'.
    kwargs:
        additional arguments to be passed to create_gaussian_input function.

    Return
    ======
    (str) name of the output file.
    """
    methods = ['scale_distance', 'increase_distance',
               'increase_distance_with_constraints']
    if method not in methods:
        raise ValueError("Non-recognized stretching method. To see the "
                         "options, check 'sith change_distance -h'")

    deltad = float(deltad)
    # Read previus file
    atoms = read(inp)
    atoms = Alignment.align_with_components(atoms)
    manipulator = MoleculeSetter(atoms)
    if deltad != 0:
        cons = np.loadtxt(file_cons, usecols=[0, 1], dtype=int)
        if type(cons[0]) is np.int64:
            cons = np.array([list(cons)])
        manipulator.xy_alignment(cons[0][0], cons[0][1], center=cons[0][0])
        eval(f'manipulator.{method}(cons, deltad)')

    manipulator.create_gaussian_input(out=out, charge=charge, **kwargs)

    return f"{out}.com"


# add2executable
def extract_bonds(readable_file):
    """
    Guesses the bonds in a molecule according to
    ase.neighborlist.natural_cutoffs.

    Parameters
    ==========
    readable_file: str
        string to the configuration file in any of the ASE readable formats.

    Return
    ======
    (list) [#bonds x 2(int)] Bonds in the molecule.

    E.g.
    sith extract_bonds optimization.xyz
    """
    atoms = read(readable_file)
    ana = Analysis(atoms)
    bonds_zm = ana.unique_bonds

    bonds = []
    for i, pairs in enumerate(bonds_zm[0]):
        if len(pairs) != 0:
            for pair in pairs:
                bonds.append([i + 1, pair + 1])
    return bonds


# add2executable
def diff_bonds(conf1, conf2, frozen_dofs='frozen_dofs.dat'):
    """
    This function returns the bonds that are in one conf but not in the other.

    Parameters
    ==========
    conf11: str
        first configuration file to be compared.
    conf2: str
        second configuration file to be compared.
    frozen_dofs: str. Default="frozen_dofs.dat"
        file with the frozen DOFs. If a rupture is obtained, it will be
        added to this file.

    Return
    ======
    (list) [#brocken_bonds x 2(int)] pair of brocken bonds

    Note
    ====
    The comparison is only in one direction, namely, this function does not
    find the bonds in conf2 that are not in conf1.
    """

    bonds1 = extract_bonds(conf1)
    bonds2 = extract_bonds(conf2)

    different_bonds = []
    for bond in bonds1:
        if bond not in bonds2:
            different_bonds.append(bond)

    with open(frozen_dofs, 'a') as fdofs:
        for bond in different_bonds:
            for index in bond:
                fdofs.write(str(index) + ' ')
            fdofs.write('F\n')
    return different_bonds


# add2executable
def conf2pdb(confile, pdbtemplate, pdboutput=None):
    """
    Transform a configuration in a file (xyz, log...) into a pdb using a pdb
    file as template.

    Parameters
    ==========
    confile: str
        path to the config file to be transformed to pdb.
    pdbtemplate: str
        path to the pdb template for the output.
    pdboutput: str. Default=None
        name of trasnformed config file with pdb format. The default name is
        the same than the confile but with pdb extension.

    Return
    ======
    (str) The name of the output pdb file.

    Note
    ====
        the pdb file must contain the same atoms in the same order than the xyz
        file, this file only would change the coordinates.

    E.g.
    sith conf2pdb optimization.log template.pdb
    sith conf2pdb optimization.xyz template.pdb
    """
    if pdboutput is None:
        pdboutput = confile.split('.')[0] + '.pdb'

    atoms_ref = read(pdbtemplate)
    atoms_xyz = read(confile)
    assert len(atoms_ref) == len(atoms_xyz), "The number of atoms in the" + \
        " reference and the template does not coincide"
    assert \
        atoms_ref.get_chemical_symbols() == atoms_xyz.get_chemical_symbols(), \
        "The atoms in the reference and the template does not coincide"
    new_positions = atoms_xyz.positions.copy()

    atoms_ref.set_positions(new_positions)
    atoms_xyz.arrays['residuenames'] = atoms_ref.arrays['residuenames']
    atoms_xyz.arrays['residuenumbers'] = atoms_ref.arrays['residuenumbers']
    atoms_xyz.arrays['atomtypes'] = atoms_ref.arrays['atomtypes']

    write(pdboutput, atoms_ref)

    return pdboutput


# add2executable
def all_xyz2pdb(template, output_patern=None, xyzdir=''):
    """
    Transform all xyz files of the directory where it is executed into a pdb
    using a pdb file as template.

    Parameters
    ==========
    pdbtemplate: str
        path to the pdb template for the output.
    output_patern: str. Default=None
        the name of the output will be <this string>-<n>.pdb, where is is an
        increasing index, from 1 to the number of xyz files.
    xyzdir: str. Default=""
        directory containing all the xyz files you want to transform.

    Return
    ======
    (str) The name of the outputs of each pdb file.

    Note
    ====
        the pdb file must contain the same atoms in the same order than the xyz
        file, this file only would change the coordinates.

    E.g.
    cd dir_with_xyz_files ; sith all_xyz2pdb
    """
    configs = glob.glob(xyzdir + '*.xyz')
    configs.sort()
    n = 1
    outfiles = []
    for config in configs:
        if output_patern is None:
            pdboutput = config.split('.')[0] + '.pdb'
        else:
            pdboutput = output_patern + f'-{n}.pdb'
            print(pdboutput)
            n += 1
        outfiles.append(conf2pdb(config, template, pdboutput=pdboutput))
    return outfiles


def all_hydrogen_atoms(mol):
    """
    Finds the indexes of all the hydrogen atoms in the peptide from Atoms
    object.

    Parameters
    ==========
    mol: string to config file or ase.Atoms object
        ASE.Atoms object to extract the hydrogen indexes.

    Return
    ======
    (list) [#h_atoms(int)] Indexes corresponding to Hydrogen atoms in 1-based
    index.
    """
    if type(mol) is str:
        mol = read(mol)
    indexes = np.where(mol.get_atomic_numbers() == 1)[0]

    return indexes + 1


# add2executable
def distance(file, index1, index2):
    """
    Computes the distance between two atoms in the last configuration
    of a trajectory file (e.g. .log file from gaussian).

    Parameters
    ==========
    file: str
        name of the file that contains the trajectory.
    index1: int
        index of the first atom to compute distances
    index2: int
        index of the second atom to compute distances

    Return
    ======
    (float)  Distance between atoms corresponding with atom with index1 and
    index2.

    E.g.
    sith distance optimization.log 1 20
    """
    index1 = int(index1)
    index2 = int(index2)
    atoms = read(file)
    d = atoms.get_distance(index1 - 1, index2 - 1)

    return d


# add2executable
def F_max_stretch(ds, pep, fd='frozen_dofs.dat'):
    """
    Compute the force in the end atoms of the configuration of
    maximum stretching.

    Parameters
    ==========
    ds: str
        path to the directory containing the peptide.
    pep: str
        name of the peptide on one letter amino acid code.
    fd: str. Default=frozen_dofs.dat
        name of the file containing the constrained atoms in the first
        line.

    Return
    ======
    (float) Force applied at the extremes of the peptide for the
    maximum stretched configuration.
    """

    fd = ds + pep + '/' + fd
    with open(fd, 'r') as files:
        lines = files.read()
        dof = int(lines.split()[0])

    # find the last stretched conf
    all_files = glob.glob(ds + pep + '/' + pep + '*.log')
    all_files.sort()
    last_conf = all_files[-1]
    atoms = read(last_conf)
    force_mag = np.linalg.norm(atoms.get_forces()[dof - 1])

    return force_mag * 960  # KJmol-1/nm


def F_stretch(logfile, index):
    """
    Extracts the magnitud of the force applied to an atoms of the configuration
    of maximum stretching.

    Parameters
    ==========
    logfile:
        file containing the optmimization or force calculation info.
    index:
        index of the atom that you want to extract the force. 1-based (like
        gaussian).

    Return
    ======
    (float) Force applied at the selected atom.
    """
    dof = index
    conf = logfile

    atoms = read(conf)
    force_mag = np.linalg.norm(atoms.get_forces()[dof - 1])

    return force_mag * 960  # KJmol-1/nm
