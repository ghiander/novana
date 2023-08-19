from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from novana.exception import NovanaError


class MoleculeBuilder(object):
    @staticmethod
    def create_rwmol_from_smiles(smiles, flatten_mol):
        """Flattens and creates an editable molecule."""
        rwmol = MoleculeBuilder._create_rwmol_from_smiles(smiles, flatten_mol)
        if rwmol is None:
            raise NovanaError("The molecule was not created for "
                              f"the input '{smiles}'")
        return rwmol

    @staticmethod
    def _create_rwmol_from_smiles(smiles, flatten_mol):
        """Business logic for rwmol creation."""
        if isinstance(smiles, str):
            if flatten_mol:
                smiles = MoleculeBuilder._remove_stereochemistry_from_smiles(
                    smiles)
            mol = Chem.MolFromSmiles(smiles)
            mol = MoleculeBuilder._keep_largest_fragment(mol)
            if mol is None:
                return None
            return Chem.RWMol(mol)

    @staticmethod
    def _remove_stereochemistry_from_smiles(smiles):
        return smiles.replace('@', '')

    @staticmethod
    def _keep_largest_fragment(mol):
        """
        Returns the largest fragment in a mixture.
        The fragment should have at least a ring,
        otherwise None will be returned.

        """
        mols = list(Chem.rdmolops.GetMolFrags(mol, asMols=True))
        if mols:
            mols.sort(reverse=True, key=lambda m: m.GetNumAtoms())
            for m in mols:
                if rdMolDescriptors.CalcNumRings(m) > 0:
                    return m
            # If no molecules with rings
            return None
