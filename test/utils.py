from rdkit import Chem


def mol_to_smiles(mol):
    """Converts an RDKit molecule to SMILES."""
    return Chem.MolToSmiles(mol)


def mol_to_inchikey(mol):
    """Converts an RDKit molecule to InChI Key."""
    return Chem.MolToInchiKey(mol)
