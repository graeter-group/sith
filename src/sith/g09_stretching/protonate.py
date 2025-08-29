# Note: This code protonates the protein to make it neutral in charge.
import numpy as np
from sith.utils.peptides import PepSetter
from ase.visualize import view
from ase import Atom
from ase.io import write


class Protonize:
    """
    Protonize charged amino acids to make them more stable in QM.

    Parameters
    ==========
    pdb: str.
        pdb file of the peptide.
    """
    def __init__(self, pdb):
        self.pepset = PepSetter(pdb)
        self.atoms = self.pepset.atoms
        all_aa = np.array(list(self.pepset.amino_name.values()))

        # K: Lysine, R: Arginine, D: Aspartic, E: Glutamine
        charged_aa = ['LYS', 'ARG', 'ASP', 'GLU']

        # identify which are the indexes of the charged amino acids
        self.matched_aa = {}
        for aa in charged_aa:
            aa_indexes = np.where(all_aa == aa)[0] + 1
            self.matched_aa[aa] = aa_indexes
        self.todel = []
        self.atoms2add = []
        self.residuenames2add = []
        self.residuenumbers2add = []
        self.atomtypes2add = []

    # Negative
    def asp(self):
        """
        In case ASP

        Return
        ======
        (ase.Atoms) atoms to be added.
        """
        for aa_index in self.matched_aa['ASP']:
            od1 = self.pepset.amino_info[aa_index]['OD1'] - 1
            cg = self.pepset.amino_info[aa_index]['CG'] - 1
            new_o = 1.8 * self.pepset.atoms[od1].position - \
                0.8 * self.pepset.atoms[cg].position
            atom = Atom('H', new_o)
            self.atoms2add.append(atom)
            self.residuenames2add.append('ASP')
            self.residuenumbers2add.append(aa_index)
            self.atomtypes2add.append('HO')
        return self.atoms2add

    def glu(self):
        """
        In case of GLU

        Return
        ======
        (ase.Atoms) atoms to be added.
        """
        for aa_index in self.matched_aa['GLU']:
            od1 = self.pepset.amino_info[aa_index]['OE1'] - 1
            cg = self.pepset.amino_info[aa_index]['CD'] - 1
            new_o = 1.8 * self.pepset.atoms[od1].position - \
                0.8 * self.pepset.atoms[cg].position
            atom = Atom('H', new_o)
            self.atoms2add.append(atom)
            self.residuenames2add.append('GLU')
            self.residuenumbers2add.append(aa_index)
            self.atomtypes2add.append('HO')
        return self.atoms2add

    # Positive
    def arg(self):
        """
        In case of ARG

        Return
        ======
        (ase.Atoms) atoms to be added.
        """
        for aa_index in self.matched_aa['ARG']:
            h = self.pepset.amino_info[aa_index]['2HH2'] - 1
            self.todel.append(h)
        return self.todel

    def lys(self):
        """
        In case LYS

        Return
        ======
        (ase.Atoms) atoms to be added.
        """
        for aa_index in self.matched_aa['LYS']:
            try:
                h = self.pepset.amino_info[aa_index]['HZ3'] - 1
            except KeyError:
                h = self.pepset.amino_info[aa_index]['3HZ'] - 1
            self.todel.append(h)
        return self.todel

    def create_atoms(self):
        """
        Add the atoms.

        Return
        ======
        (None)
        """
        # Populate Atoms to add
        self.asp()
        self.glu()

        # Populate "to delete"
        self.arg()
        self.lys()

        del self.atoms[self.todel]

        for i, Hat in enumerate(self.atoms2add):
            self.atoms += Hat
            self.atoms.arrays['residuenames'][-1] = self.residuenames2add[i]
            self.atoms.arrays['atomtypes'][-1] = self.atomtypes2add[i]
            self.atoms.arrays['residuenumbers'][-1] = \
                self.residuenumbers2add[i]


# add2executable
def protonate(pdb, output):
    """
    Add or remove the H atoms neccessary to neutralize charges in the atoms
    charged amino acids.

    Parameters
    ==========
    pdb: str
        name of the pdb file with the amino acids.
    output:
        name of the output file.

    Return
    ======
    (ase.Atoms) new structures.
    """
    prot = Protonize(pdb)
    prot.create_atoms()
    write(output, prot.atoms)

    return prot.atoms
