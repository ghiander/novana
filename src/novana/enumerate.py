from novana.core import adjust_hydrogens
from novana.core import adjust_atom_hs_before_split
from novana.process import TerminalAtomRemover
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import rdmolops


def find_removable_groups(rwmol):
    initial_groups = _find_terminal_atom_groups(rwmol)
    expanded_groups = _expand_atom_groups(rwmol, initial_groups)
    return deduplicate_unsorted_list_of_lists(expanded_groups)


def _find_terminal_atom_groups(rwmol):
    """Create a list of singletons of terminal atoms in the molecule."""
    terminal_atoms_groups = []
    for bond in list(rwmol.GetBonds()):
        for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
            if TerminalAtomRemover._is_atom_single_bonded(atom):
                terminal_atoms_groups.append([atom.GetIdx()])
    return terminal_atoms_groups


def _expand_atom_groups(rwmol, atoms_groups):
    for group in atoms_groups:
        for atom_idx in group:
            atom = rwmol.GetAtomWithIdx(atom_idx)
            bonds = atom.GetBonds()
            for b in bonds:
                neighbour_atom = _get_neighbour_atom(b, atom_idx)
                if neighbour_atom.IsInRing():
                    continue
                neighbour_atom_idx = neighbour_atom.GetIdx()
                if neighbour_atom_idx not in group:
                    group.append(neighbour_atom_idx)
    return atoms_groups


def _get_neighbour_atom(bond, query_atom_idx):
    if bond.GetBeginAtom().GetIdx() == query_atom_idx:
        return bond.GetEndAtom()
    else:
        return bond.GetBeginAtom()


def deduplicate_unsorted_list_of_lists(list_of_lists):
    return [list(x) for x in set(tuple(sorted(lst)) for lst in list_of_lists)]


def remove_isolated_atoms(rwmol):
    isolated_atoms = [atom.GetIdx() for atom in rwmol.GetAtoms() if atom.GetDegree() == 0]
    for idx in sorted(isolated_atoms, reverse=True):
        rwmol.RemoveAtom(idx)
    return rwmol


def enumerate_substructures(rwmol, removable_groups):
    """Enumerates the partial decompositions of a molecule."""
    enumerated_mols = list()
    for group in removable_groups:
        rw_substructure = create_molecule_with_removed_group(rwmol, group)
        # If decomposition has created multiple fragments, skip it
        if not _has_molecule_multiple_fragments(rw_substructure):
            enumerated_mols.append(rw_substructure)
    return enumerated_mols


def create_molecule_with_removed_group(rwmol, group):
    """Creates a new molecule with the specified group removed."""
    _rwmol = deepcopy(rwmol)
    _rwmol = _delete_group_from_mol(_rwmol, group)
    _rwmol = adjust_hydrogens(_rwmol)
    Chem.SanitizeMol(_rwmol)
    return _rwmol


def _delete_group_from_mol(rwmol, group):
    """Delete a group of atoms and their bonds from a molecule."""
    reversed_group = sorted(group, reverse=True)  # To avoid index shifting
    for atom_idx in reversed_group:
        bonds = rwmol.GetAtomWithIdx(atom_idx).GetBonds()
        for bond in bonds:
            neighbour_atom = _get_neighbour_atom(bond, atom_idx)
            adjust_atom_hs_before_split(neighbour_atom, bond)
            rwmol.RemoveBond(atom_idx,
                             neighbour_atom.GetIdx())
    return remove_isolated_atoms(rwmol)


def _has_molecule_multiple_fragments(mol):
    """Check if a molecule has multiple fragments."""
    fragment_ids = rdmolops.GetMolFrags(mol)
    return len(fragment_ids) > 1
