import numpy as np
from ase.io import read
import sys
import matplotlib.pyplot as plt
from sith.utils.miscellaneous import output_terminal
from ase.visualize import view
from sith.utils.molecules import MoleculeSetter, Alignment
from sith.utils.peptides import PepSetter
from ase import Atoms
from ase.io import write
import glob


# add2executable
def info_from_opt(logfile, pdb_reference, pattern):
    """
    Extracts configurations from a gaussian optimization trajectory and removes
    those outliers in the total energy.

    Parameters
    ==========
    logfile: str
        gaussian log file.
    pdb_reference: str
        pdb with the information of the peptide.
    pattern: str
        pattern of the output. The results are xyz files called
        <pattern><i>.xyz

    Return
    ======
    (ase.Atoms) new trajectory without the outliers. it creates the xyz files
    in between.
    """
    # read configurations
    atoms = read(logfile, index=':')

    # read energies
    energies = output_terminal("grep 'SCF Done:' " + logfile
                               + " | awk '{print $5}'",
                               print_output=False)
    energies = np.array(energies.split('\n')[:-1], dtype=float)

    # read amino acids info, recognize extreme indexes
    ps = PepSetter(pdb_reference)
    res = list(ps.amino_info.keys())
    middle = len(res) // 2 + 1
    ind1 = ps.amino_info[res[0]]['CH3'] - 1
    ind2 = ps.amino_info[res[-1]]['CH3'] - 1

    # choose third atom for orientation.
    if ps.amino_name[middle] != 'GLY':  # Glycine does not have CB
        ind3 = ps.amino_info[3]['CB'] - 1
    else:  # In case of 3 Glycine
        ind3 = ps.amino_info[middle]['CA'] - 1

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

    # align all the structures in the same plane. This guarantee that
    # the average of two structures is the intermedia structure between them
    all_atoms = [Alignment.align_with_components(conf) for conf in atoms]
    for conf in all_atoms:
        ms = MoleculeSetter(conf)
        ms.xy_alignment(ind1, ind2, ind3)

    # now make the trajectory continuos. if the distance between extremes is
    # larger than 0.2A, configurations with the intermedia distances are
    # created.
    new_set = [atoms[0]]
    i = 1
    while i < len(atoms):
        conf = atoms[i]
        di = new_set[-1].get_distance(ind1, ind2)
        df = conf.get_distance(ind1, ind2)
        deltad = df - di
        if abs(deltad) > 0.2:
            n_intermedia = int(abs(deltad / 0.2))
            fraction = 1/(n_intermedia + 1)
            inbetween = []
            for n in range(1, n_intermedia + 1):
                at = Atoms(conf.get_chemical_symbols(),
                           positions = (1 - fraction * n) *
                                       new_set[-1].positions +
                                       (fraction * n) * conf.positions)
                inbetween.append(at)
            new_set.extend(inbetween)
        new_set.append(conf)
        i += 1

    # in case of last configurations does not belong to new_set, it is added
    if np.any(atoms[-1].positions != new_set[-1].positions):
        new_set.append(atoms[-1])

    # write all the trajectory
    for i, atoms in enumerate(all_atoms):
        write('{}{:03d}.xyz'.format(pattern, i), atoms)

    return all_atoms


# add2executable
def reduce_structs(dir, pattern):
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

    # jump to furthest repeated structure and save the new order in new_set
    new_set = []  # it will contain the indices of the selected indexes
    j = 0
    while j < len(all_files) - 1:
        new_set.append(j)
        d_ij = (all_dofs[j] - all_dofs[j + 1:])
        # rescale distances and angles to have them in the same order of
        # magnitude. I could also evaluate the approximation for distances and
        # then for angles and use logicaland.
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

    # make it continuous
    subdofs = all_dofs[new_set]
    search_intermedia = True
    while search_intermedia:
        search_intermedia = False
        # make distances continuous
        d_ij = subdofs[1:] - subdofs[:-1]
        d_ij = delta_angles_continuous(d_ij, nrs)
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
    output_terminal("if [ ! -d subset ]; then mkdir subset; fi")
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
