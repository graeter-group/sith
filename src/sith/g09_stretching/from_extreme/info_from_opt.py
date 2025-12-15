import numpy as np
from ase.io import read
from sith.utils.miscellaneous import output_terminal
from sith.utils.molecules import Alignment, MoleculeSetter
from ase.io import write
from ase import Atoms
import os
import glob


# add2executable
def info_from_opt(logfile, pattern):
    """
    Extracts configurations from a gaussian optimization trajectory and removes
    the structures that are outliers in the total energy.

    Parameters
    ==========
    logfile: str
        gaussian log file.
    pattern: str
        pattern of the output. The results are xyz files called
        <pattern><i>.xyz

    Return
    ======
    (ase.Atoms) new trajectory without the outliers. It creates the xyz files
    in between.

    Note
    ====
    You can use
    :func:`sith.g09_stretching.from_extreme.info_from_opt.reduce_struct` to
    make the DOFs of output trajectory continuous.
    """
    # read configurations
    atoms = read(logfile, index=':')

    try:
        energies = [conf.get_potential_energy() for conf in atoms]
        energies = np.array(energies)
    except:
        # read energies
        energies = output_terminal("grep 'SCF Done:' " + logfile
                                + " | awk '{print $5}'",
                                print_output=False)
        energies = np.array(energies.split('\n')[:-1], dtype=float)

    # remove configurations that goes up in energy. keep those that goes
    # down only. assuming local optimization
    i = 0
    while True:
        de = energies[1:] - energies[:-1]
        toremove = np.where(de > 0)[0]
        if len(toremove) == 0:
            break

        for index in toremove[::-1]:
            atoms.pop(index + 1)
        energies = np.delete(energies, toremove + 1)
        i += 1

    # align all the structures in the same plane. This guarantees that
    # the average of two structures is the intermedia structure between them
    all_atoms = [Alignment.align_with_components(conf) for conf in atoms]

    # write all the trajectory
    for i, atoms in enumerate(all_atoms[::-1]):
        write('{}{:03d}.xyz'.format(pattern, i), atoms)

    return all_atoms


# add2executable
def continuous_e2e(pattern, index1, index2, dir_output='e2e_continuous'):
    """
    From a set of configurations in xyz files with pattern
    <all><pattern><all>.xyz, it guarantees that consecutive changes in the
    distance between the atoms with index1 and index2 are lower than 0.2A.

    Parameters
    ==========
    pattern: str
        pattern of the input files. All the files with <all><pattern><all>.xyz
        will be read.
    index1: int
        index of the first atom to align in the x axis. 1-base indexing.
    index2: int
        index of the second atom to align in the x axis. 1-base indexing.
    dir_output: str. Default='e2e_continuous'
        directory where to save the output files.

    Return
    ======
    (ase.Atoms) new trajectory with continuous distance between index1 and
    index2. It creates the xyz files in dir_output.
    """
    all_atoms = []
    files = glob("*{pattern}*.xyz")
    files.sort()
    for file in files:
        atoms = read(file)
        all_atoms.append(atoms)

    # choose randomlly an atom betweeen index1 and index2 that is not H and not
    # in the x axis
    # ==== align all the structures with index1 and index2 in the x axis.
    for i, conf in enumerate(all_atoms):
        ms = MoleculeSetter(conf)
        ms.xy_alignment(index1, index2)
        all_atoms[i] = ms.atoms
    # ==== inbetween
    positions = np.array([conf.get_positions() for conf in all_atoms])
    inbetween = np.logical_and(positions[:, :, 0].T
                               > positions[:, index1 - 1, 0],
                               positions[:, :, 0].T
                               < positions[:, index2 - 1, 0])
    inbetween = inbetween.T
    inbetween = np.all(inbetween, axis=0)
    # ==== avoid
    d_xaxis = np.sqrt(positions[:, :, 0]**2 + positions[:, :, 1]**2)
    not_xaxis = np.all(d_xaxis != 0, axis=0)  # A
    # ==== heavy atoms
    heavys = np.array(all_atoms[0].get_chemical_symbols()) != 'H'
    # ==== random choice
    condition = np.logical_and(heavys, not_xaxis)
    condition = np.logical_and(condition, inbetween)
    subset = np.where(condition)[0]
    if len(subset) == 0:
        raise ValueError("No heavy atom found between index1 and index2 that"
                         "is not in the x axis")
    index3 = np.random.choice(subset) + 1

    # alignment with 3 atoms. All the configurations are aligned with the
    # same 3 atoms. This guarantees that the average of two structures is the
    # intermedia structure between them
    for i, conf in enumerate(all_atoms):
        ms = MoleculeSetter(conf)
        ms.xy_alignment(index1, index2, index3)
        all_atoms[i] = ms.atoms

    # now make the trajectory continuos. if the distance between extremes is
    # larger than 0.2A, configurations with the intermedia distances are
    # created.
    new_set = [all_atoms[0]]
    i = 1
    while i < len(all_atoms):
        conf = all_atoms[i]
        di = new_set[-1].get_distance(index1 - 1, index2 - 1)
        df = conf.get_distance(index1 - 1, index2 - 1)
        deltad = df - di
        if abs(deltad) > 0.2:
            n_intermedia = int(abs(deltad / 0.2))
            fraction = 1 / (n_intermedia + 1)
            inbetween = []
            for n in range(1, n_intermedia + 1):
                at = Atoms(conf.get_chemical_symbols(),
                           positions=((1 - fraction * n)
                                      * new_set[-1].positions
                                      + (fraction * n) * conf.positions))
                inbetween.append(at)
            new_set.extend(inbetween)
        new_set.append(conf)
        i += 1

    # in case of last configurations does not belong to new_set, it is added
    if np.any(all_atoms[-1].positions != new_set[-1].positions):
        new_set.append(all_atoms[-1])

    # write all the trajectory
    os.makedirs(dir_output, exist_ok=True)
    for i, atoms in enumerate(new_set):
        write("{}/{}{:03d}.xyz".format(dir_output, pattern, i), atoms)

    return all_atoms


# add2executable
def reduce_structs(dir, pattern, print_reduced=False):
    """
    Check all the <all>-dofs.dat files and remove those files that represent
    irrelevant changes. It creates intermedias to guarantee continuous dofs.

    Parameters
    ==========
    dir: str
        path to the directory containing all the <all>-dofs.dat files.
    pattern: str
        pattern for the output files, all of them will be then
        <pattern><i>.dat

    Return
    ======
    (numpy.array) new set of DOFs [structure][DOF], which are reduced. In
    between, it also creates all <pattern><i>,dat.
    """
    # find dofs files
    all_files = glob.glob(f"{dir}/*-dofs.dat")
    all_files.sort()

    # extract dofs definitions and dofs values
    dofs_ref = np.loadtxt(all_files[0], delimiter='=',
                          comments='      Variables:', usecols=0, dtype=str)
    nrs = len([r for r in dofs_ref if r[1] == 'R'])
    all_dofs = []
    for file in all_files:
        dofs = np.loadtxt(file, delimiter='=',
                          comments='      Variables:', usecols=0, dtype=str)
        assert (dofs == dofs_ref).all(), \
            f"{file} has different dofs than {all_files[0]}"
        dofs = np.loadtxt(file, delimiter='=',
                          comments='      Variables:', usecols=1)
        all_dofs.append(dofs)
    all_dofs = np.array(all_dofs)  # [struct][dof]

    # jump to the furthest repeated structure and save the new order in new_set
    new_set = []  # it will contain the indices of the selected indexes
    j = 0
    while j < len(all_files) - 1:
        new_set.append(j)
        d_ij = (all_dofs[j] - all_dofs[j + 1:])
        # rescale distances to have them in the same order of magnitude of the
        # angles.
        d_ij = delta_angles_continuous(d_ij, nrs)
        d_ij[:, nrs:] *= 1e-3  # trans 1 degree
        d_ij = abs(d_ij)
        close = np.isclose(d_ij, 0, atol=1e-3)
        by_struc = np.all(close, axis=1)
        same = np.where(by_struc)[0] + j
        if len(same) != 0:
            j = same[-1]
        j += 1
    if new_set[-1] != len(all_files) - 1:
        new_set.append(len(all_files) - 1)
    
    if print_reduced:
        for i in new_set:
            print(all_files[i])

    # make it continuous
    subdofs = all_dofs[new_set]
    search_intermedia = True
    while search_intermedia:
        search_intermedia = False
        # make distances continuous
        d_ij = subdofs[1:] - subdofs[:-1]
        d_ij = delta_angles_continuous(d_ij, nrs)
        # Create intermedias between distances larger than 0.05A or angles
        # larger than 10 degrees.
        condition1 = np.any(abs(d_ij[:, :nrs]) > 0.05, axis=1)
        condition2 = np.any(abs(d_ij[:, nrs:]) > 10, axis=1)
        toModify = np.where(np.logical_or(condition1,
                                          condition2))[0]
        for i in toModify[::-1]:
            search_intermedia = True
            # make angles continuous
            subdofs = np.insert(subdofs, i + 1,
                                d_ij[i] / 2 + subdofs[i],
                                axis=0)

    # copy relevant files to a directory called subset
    output_terminal("mkdir -p subset")
    for i, struct in enumerate(subdofs):
        with open('./subset/{}{:03d}.dat'.format(pattern,
                                                 i), "w") as i_struct_file:
            for j, dof in enumerate(struct):
                i_struct_file.write(f'{dofs_ref[j]}={dof}\n')

    return subdofs


def delta_angles_continuous(d_ij, nrs):
    """
    Make the difference of the angles continuous considering the periodicity of
    them.

    Parameters
    ==========
    d_ij: numpy.Array
        set of delta dofs with form [structure][DOF]
    nrs: int
        Number of distances.

    Return
    ======
    (numpy.Array) new set of delta DOFs with delta angles continuous.

    Note
    ----
    Here, it is assumed that the first <nrs> DOFs in d_ij are distances and
    the rest are angles.
    """
    condition = d_ij[:, nrs:] < -180
    while condition.any():
        d_ij[:, nrs:][condition] += 360
        condition = d_ij[:, nrs:] < -180
    condition = d_ij[:, nrs:] > 180
    while condition.any():
        d_ij[:, nrs:][condition] -= 360
        condition = d_ij[:, nrs:] > 180
    return d_ij
