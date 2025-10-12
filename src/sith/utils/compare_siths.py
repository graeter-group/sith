import numpy as np


def index_dof(sith, target: tuple):
    """
    Search the index of a specific dof.

    Parameters
    ==========
    sith: sith.SITH
        SITH object.
    target: np.ndarray
        degree of freedom.

    Return
    ======
    (int) index
    """
    for i, dof in enumerate(sith.dim_indices):
        dof_wo_0 = dof[np.nonzero(dof)[0]]
        target_wo_0 = target[np.nonzero(target)[0]]

        if len(dof_wo_0) != len(target_wo_0):
            continue

        if (dof_wo_0 == target_wo_0).all() or \
           (dof_wo_0 == target_wo_0[::-1]).all():
            return i
    raise ValueError("Non-found dof.")


def organize_lengths(sith1, sith2):
    """
    sort the distances in sith1 such that they fit with sith2.
    # TODO: confirm that the description is correct and it is
    # not the other way arround

    Parameters
    ==========
    sith1: sith.SITH
        obtject of reference.
    sith2: sith.SITH
        object to be compared.

    Return
    ======
    (numpy.array) set of sorted indexes that makes the distance of sith2 to
    coincide with the distances of sith1.
    """
    n_distances = sith2.dims[1]

    reorder = []
    for dof in np.array(sith1.dim_indices[:n_distances]):
        new_index = index_dof(sith2, dof)
        reorder.append(new_index)

    return np.array(reorder)


def extract_common(sith1, sith2):
    """
    Finds the DOfs that are in common in both sith objects

    Parameters
    ==========
    sith1: sith.SITH
        obtject of reference.
    sith2: sith.SITH
        object to be compared.

    Return
    ======
    (tuple), list of indexes in sith1 that are in sith2, organized list.
    """
    sith1_indx = []
    sith2_indx = []
    for i, dof in enumerate(np.array(sith1.dim_indices)):
        try:
            new_index = index_dof(sith2, dof)
            sith1_indx.append(i)
            sith2_indx.append(new_index)
        except ValueError:
            continue

    return np.array(sith1_indx), np.array(sith2_indx)


def extract_diff(sith1, sith2):
    """
    Finds the indexes of the DOfs that are in sith1 that are not in sith2 and
    the other way around.

    Parameters
    ==========
    sith1: sith.SITH
        object to compare.
    sith2: sith.SITH
        object to compare.

    Return
    ======
    (tuple), list of indexes in sith1 that are in sith2, organized list.
    """
    sith1_indx = []
    sith2_indx = []
    for i, dof in enumerate(np.array(sith1.dim_indices)):
        try:
            index_dof(sith2, dof)
        except ValueError:
            sith1_indx.append(i)

    for i, dof in enumerate(np.array(sith2.dim_indices)):
        try:
            index_dof(sith1, dof)
        except ValueError:
            sith2_indx.append(i)
            continue

    return np.array(sith1_indx), np.array(sith2_indx)
