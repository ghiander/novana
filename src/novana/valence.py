from math import ceil
from rdkit import Chem
from novana.utils import get_logger
from novana.utils import disable_rdkit_logging


# Configure logging
log = get_logger(__name__)
disable_rdkit_logging()

# Static objects
periodic_table = Chem.GetPeriodicTable()


class AtomProperties(object):
    def __init__(self, atom):
        self.atom = atom
        self.symbol = atom.GetSymbol()
        self.index = atom.GetIdx()
        self.atomic_num = atom.GetAtomicNum()
        self.explict_hs_num = atom.GetNumExplicitHs()
        self.total_valence = atom.GetTotalValence()
        self.formal_charge = atom.GetFormalCharge()
        self.default_valence = self._calculate_default_valence(periodic_table)

    def _calculate_default_valence(self, periodic_table):
        return periodic_table.GetDefaultValence(self.atomic_num)

    def are_there_explicit_hs(self):
        if self.explict_hs_num > 0:
            return True
        return False

    def log_debug_properties(self):
        log.debug(f"Atom: {self.symbol}")
        log.debug(f"Atom index: {self.index}")
        log.debug(f"Charge: {self.formal_charge}")
        log.debug(f"Total valence: {self.total_valence}")
        log.debug(f"Default valence: {self.default_valence}")
        log.debug(f"Explicit hydrogens {self.explict_hs_num}")


class HydrogenAdjuster(object):
    def __init__(self, mol):
        self.mol = mol

    def process_and_return(self):
        self.run()
        return self.mol

    def run(self):
        """
        Sanitisation heuristics to fix hypervalent atoms
        by removing hydrogens in excess; and hypovalent
        atoms by adding hydrogens. Both processes are done
        in order to achieve the default valence of the atom.

        """
        log.debug(Chem.MolToSmiles(self.mol))
        for bond in list(self.mol.GetBonds()):
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                atom_props = AtomProperties(atom)
                if atom_props.are_there_explicit_hs:
                    bond_type = bond.GetBondTypeAsDouble()
                    atom_props.log_debug_properties()
                    log.debug(f"Bond type as double: {bond_type}")

                    # If total valence is greater than default, hydrogens
                    # need to be removed and charges neutralised
                    if atom_props.total_valence > atom_props.default_valence:
                        self.remove_hydrogens(self, atom, atom_props,
                                              bond_type)

                    # If total valence is smaller than default, hydrogens
                    # need to be added and charges neutralised
                    # In addition, if aromatic bonds are present,
                    # they should be converted into single to accommodate
                    # the previous changes, otherwise kekulisation will fail
                    elif atom_props.total_valence < atom_props.default_valence:
                        self.add_hydrogens(atom, atom_props, bond, bond_type)

                    atom.UpdatePropertyCache()
                    break

    @staticmethod
    def remove_hydrogens(self, atom, atom_props, bond_type):
        hs_to_remove = int(
            ceil((atom_props.total_valence - atom_props.default_valence)
                 / bond_type))
        if hs_to_remove <= atom_props.explict_hs_num and hs_to_remove > 0:
            atom.SetNumExplicitHs(atom_props.explict_hs_num - hs_to_remove)
            atom.SetFormalCharge(0)
        # TODO: Restore
        # Chem.SanitizeMol(self.mol)
        # Chem.rdmolops.SanitizeFlags.SANITIZE_NONE

    def add_hydrogens(self, atom, atom_props, bond, bond_type):
        """"""
        if bond_type == 1.5:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
        atom.SetNumExplicitHs(atom_props.explict_hs_num + 1)
        atom.SetFormalCharge(0)
        Chem.SanitizeMol(self.mol)
        Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
