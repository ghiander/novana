from rdkit import Chem


class TerminalAtomRemover(object):
    def __init__(self, rwmol):
        self.mol = rwmol

    def process_and_return(self):
        self.remove_atoms_recursively()
        return self.mol

    def remove_atoms_recursively(self):
        """
        Removes recursively all terminal
        atoms and their bonds.

        """
        init_atom_length = len(self.mol.GetAtoms())
        final_atom_length = 0   # dummy length to start the while
        while init_atom_length != final_atom_length:
            init_atom_length = len(self.mol.GetAtoms())
            self._remove_atoms()
            final_atom_length = len(self.mol.GetAtoms())
        Chem.SanitizeMol(self.mol)
        Chem.rdmolops.SanitizeFlags.SANITIZE_NONE

    def _remove_atoms(self):
        """
        Iterates through the molecule's atoms and
        removes all atoms that are terminal (and
        their bonds).
        """
        for bond in list(self.mol.GetBonds()):
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                if TerminalAtomRemover._is_atom_single_bonded(atom):
                    nbr_atom = bond.GetOtherAtom(atom)
                    if TerminalAtomRemover._is_atom_multi_bonded(nbr_atom):
                        self._remove_single_bonded_atom_from_mol(
                            single_bonded_atom=atom,
                            neighbour_atom=nbr_atom,
                            bond=bond)
                        break

    @staticmethod
    def _is_atom_single_bonded(atom):
        if atom.GetDegree() == 1:
            return True
        return False

    @staticmethod
    def _is_atom_multi_bonded(atom):
        if atom.GetDegree() > 1:
            return True
        return False

    def _remove_single_bonded_atom_from_mol(self, single_bonded_atom,
                                            neighbour_atom, bond):
        neighbour_atom = \
            TerminalAtomRemover._adjust_atom_hs_before_split(
                neighbour_atom, bond)
        self.mol.RemoveBond(single_bonded_atom.GetIdx(),
                            neighbour_atom.GetIdx())
        self.mol.RemoveAtom(single_bonded_atom.GetIdx())

    @staticmethod
    def _adjust_atom_hs_before_split(atom, bond):
        """
        Increase the explicit hydrogens of an
        atom by the rounded value of its bond type.

        """
        atom.SetNumExplicitHs(
            atom.GetNumExplicitHs()
            + int(bond.GetBondTypeAsDouble()))
        return atom


class SingleBondedCarborConverter(object):
    def __init__(self, rwmol):
        self.mol = rwmol

    def process_and_return(self):
        self.convert()
        return self.mol

    def convert(self):
        """Converts all atoms into single-bonded carbons."""
        for bond in list(self.mol.GetBonds()):
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                SingleBondedCarborConverter._convert_atom(atom)
                break

    @staticmethod
    def _convert_atom(atom):
        atom.SetIsAromatic(False)
        atom.SetAtomicNum(6)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
