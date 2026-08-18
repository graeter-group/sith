# Code modified from:
# https://github.com/chunxiangzheng/gaussian_log_file_converter/
import re
from ase import Atoms
from ase.io import write, read
from sith.utils.molecules import MoleculeSetter


code = {"1": "H", "2": "He", "3": "Li", "4": "Be", "5": "B",
        "6": "C", "7": "N", "8": "O", "9": "F", "10": "Ne",
        "11": "Na", "12": "Mg", "13": "Al", "14": "Si", "15": "P",
        "16": "S", "17": "Cl", "18": "Ar", "19": "K", "20": "Ca",
        "21": "Sc", "22": "Ti", "23": "V", "24": "Cr", "25": "Mn",
        "26": "Fe", "27": "Co", "28": "Ni", "29": "Cu", "30": "Zn",
        "31": "Ga", "32": "Ge", "33": "As", "34": "Se", "35": "Br",
        "36": "Kr", "37": "Rb", "38": "Sr", "39": "Y", "40": "Zr",
        "41": "Nb", "42": "Mo", "43": "Tc", "44": "Ru", "45": "Rh",
        "46": "Pd", "47": "Ag", "48": "Cd", "49": "In", "50": "Sn",
        "51": "Sb", "52": "Te", "53": "I", "54": "Xe", "55": "Cs",
        "56": "Ba", "57": "La", "58": "Ce", "59": "Pr", "60": "Nd",
        "61": "Pm", "62": "Sm", "63": "Eu", "64": "Gd", "65": "Tb",
        "66": "Dy", "67": "Ho", "68": "Er", "69": "Tm", "70": "Yb",
        "71": "Lu", "72": "Hf", "73": "Ta", "74": "W", "75": "Re",
        "76": "Os", "77": "Ir", "78": "Pt", "79": "Au", "80": "Hg",
        "81": "Tl", "82": "Pb", "83": "Bi", "84": "Po", "85": "At",
        "86": "Rn", "87": "Fr", "88": "Ra", "89": "Ac", "90": "Th",
        "91": "Pa", "92": "U", "93": "Np", "94": "Pu", "95": "Am",
        "96": "Cm", "97": "Bk", "98": "Cf", "99": "Es", "100": "Fm",
        "101": "Md", "102": "No", "103": "Lr", "104": "Rf", "105": "Db",
        "106": "Sg", "107": "Bh", "108": "Hs", "109": "Mt", "110": "Ds",
        "111": "Rg", "112": "Uub", "113": "Uut", "114": "Uuq", "115": "Uup",
        "116": "Uuh", "117": "Uus", "118": "Uuo"}


def _getEnergy(structure):
    """
    Extract the energy of the logfile.

    Parameters
    ==========
    structure: str
        block of text (lines) corresponding to a single structure from a
        gaussian log file, as accumulated by log2xyz.

    Return
    ======
    (float) DFT energy found by gaussian. Otherwise, it returns 1000.0
    """
    for line in structure.split("\n"):
        if line.startswith(" SCF Done:"):
            arr = line.split("=")
            return float(re.split(" +", arr[1].strip())[0])
    return 1000.0


def _findInList(dataList, target):
    """
    Find the first element of a list that contains a given substring.

    Parameters
    ==========
    dataList: list of str
        List of strings to search through.
    target: str
        Substring to look for in each element of dataList.

    Return
    ======
    (int) index of the first element in dataList that contains target.
    If no element contains target, -1 is returned.
    """
    for i in range(0, len(dataList)):
        if dataList[i].find(target) != -1:
            return i
    return -1


def _getCoordinates(dataList):
    """
    Extract the atomic coordinate lines from the "Standard orientation"
    block of a gaussian structure.

    Parameters
    ==========
    dataList: list of str
        Lines of a single structure block from a gaussian log file (as
        produced by splitting the text returned by _getEnergy's input on
        newlines).

    Return
    ======
    (list of str) the lines listing the atoms and their coordinates found
    right after the "Standard orientation" header (skipping the 5 header
    lines) and before the following line of dashes.
    """
    start = _findInList(dataList, "Standard orientation")
    dataList = dataList[start + 5:]
    dataList = dataList[: _findInList(dataList, "-----")]
    return dataList


# add2executable
def log2xyz(finput, foutput=None, indexes=None):
    """
    Extract the configuration of minumum energy from a .log gaussian file of an
    optimization process in a xyz file.

    Parameters
    ==========
    finput: str
        path to the log file.
    foutput: str. Default=None
        name of the output file without extension.

    Return
    ======
    (ase.Atoms) the extracted (and, if indexes is given, xy-aligned)
    optimized structure. The same structure is also written to the
    .xyz output file.

    Note: if foutput is not given, the name output will be the same than the
    input but with xyz extension.
    """
    infoBlock = ""
    optimized = False
    optimized_structure = ""

    with open(finput, "r") as fin:
        isStructure = True
        isInfo = True
        structures = []
        currentStructure = ""
        for line in fin:
            if line.startswith(" GradGrad"):
                if isInfo:
                    isInfo = False
                if currentStructure != "":
                    structures.append((_getEnergy(currentStructure),
                                       currentStructure))
                    currentStructure = ""
                isStructure = not isStructure
            elif isInfo:
                infoBlock += line
            elif isStructure:
                currentStructure += line
            else:
                if line.find("Optimized") != -1 and not optimized:
                    optimized = True
                    optimized_structure = structures[-1][1]
                    break

        if not optimized:
            if currentStructure != "":
                structures.append((_getEnergy(currentStructure),
                                   currentStructure))
            structures = sorted(structures, key=lambda item: item[0])
            optimized_structure = structures[0][1]

    dataList = optimized_structure.split("\n")
    atoms = _getCoordinates(dataList)
    mol_pos = []
    mol_anames = []
    for atom in atoms:
        arr = atom.split()
        symbol = code.get(arr[1], 'X')
        mol_pos.append([float(arr[3]), float(arr[4]), float(arr[5])])
        mol_anames.append(symbol)

    atoms = Atoms(mol_anames, mol_pos)
    ms = MoleculeSetter(atoms)
    if indexes is not None:
        if len(indexes) == 2:
            indexes.append(None)
        ms.xy_alignment(indexes[0], indexes[1], index3=indexes[2])
    atoms = ms.atoms

    if foutput:
        prefix = foutput
    else:
        prefix = finput[:-4]

    foutput = prefix + ".xyz"
    write(foutput, atoms)

    return atoms

# add2executable


def log2xyz2(finput, foutput=None, indexes=None, frame=-1):
    """
    Extract a single frame from a gaussian log file (or any file format
    supported by ase.io.read) and write it as a .xyz file.

    Parameters
    ==========
    finput: str
        path to the input file.
    foutput: str. Default=None
        name of the output file without extension. If not given, the
        output name will be the same as the input but with xyz extension.
    indexes: list. Default=None
        indexes of up to 3 atoms used to align the structure on the xy
        plane via MoleculeSetter.xy_alignment. If None, no alignment is
        performed.
    frame: int. Default=-1
        index of the frame/structure to read from finput.

    Return
    ======
    (ase.Atoms) the extracted (and, if indexes is given, xy-aligned)
    structure. The same structure is also written to the .xyz output file.
    """
    atoms = read(finput, index=frame)

    ms = MoleculeSetter(atoms)
    if indexes is not None:
        if len(indexes) == 2:
            indexes.append(None)
        ms.xy_alignment(indexes[0], indexes[1], index3=indexes[2])
    atoms = ms.atoms

    if foutput:
        prefix = foutput
    else:
        prefix = finput[:-4]

    foutput = prefix + ".xyz"
    write(foutput, atoms, format='xyz')

    return atoms
