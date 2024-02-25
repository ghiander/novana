from novana.core import adjust_hydrogens
from novana.core import convert_into_single_bonded_carbons
from novana.core import convert_into_single_bonded_atoms
from novana.core import convert_heteroatoms_to_stars
from novana.core import create_rwmol_from_smiles
from novana.core import recursively_remove_single_bonded_atoms


def scaffold_from_smiles(smiles, flatten_mol=True,
                         generalise_heteroatoms=False):
    """Returns the scaffold of a molecule."""
    rwmol = create_rwmol_from_smiles(smiles, flatten_mol)
    recursively_remove_single_bonded_atoms(rwmol)
    if generalise_heteroatoms:
        convert_heteroatoms_to_stars(rwmol)
    adjust_hydrogens(rwmol)
    return rwmol.GetMol()


def shape_from_smiles(smiles, flatten_mol=True,
                      retain_atom_types=False):
    """Returns the shape of a molecule."""
    rwmol = create_rwmol_from_smiles(smiles, flatten_mol)
    recursively_remove_single_bonded_atoms(rwmol)
    if retain_atom_types:
        convert_into_single_bonded_atoms(rwmol)
    else:
        convert_into_single_bonded_carbons(rwmol)
    return rwmol.GetMol()


def molecule_scaffold_shape_from_smiles(smiles,
                                        flatten_mol=True):
    """Returns molecule, scaffold, and shape of a molecule."""
    rwmol = create_rwmol_from_smiles(smiles, flatten_mol)
    mol = rwmol.GetMol()
    # Get scaffold and clean it
    recursively_remove_single_bonded_atoms(rwmol)
    adjust_hydrogens(rwmol)
    scaffold = rwmol.GetMol()

    # Convert scaffold into shape
    convert_into_single_bonded_carbons(rwmol)
    shape = rwmol.GetMol()
    return mol, scaffold, shape
