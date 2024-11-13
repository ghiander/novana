from rdkit import Chem


def mol_to_smiles(mol):
    """Converts an RDKit molecule to SMILES."""
    return Chem.MolToSmiles(mol)


def mol_to_inchikey(mol):
    """Converts an RDKit molecule to InChI Key."""
    return Chem.MolToInchiKey(mol)


def mol_list_to_smiles_list(mols):
    """Converts a list of RDKit molecules to SMILES."""
    return [mol_to_smiles(mol) for mol in mols]
