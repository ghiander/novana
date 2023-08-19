import pytest
from novana.api import scaffold_from_smiles
from novana.api import shape_from_smiles
from novana.exception import NovanaError


def test_exception_mixture_no_rings_scaffold():
    """Raises an error if no molecules with rings."""
    smiles = "OC(CC(=O)O)(CC(=O)O)C(=O)O.OC(CC(=O)O)(CC(=O)O)C(=O)O"
    with pytest.raises(NovanaError):
        scaffold_from_smiles(smiles)


def test_exception_molecule_no_rings_scaffold():
    """Raises an error if no molecules with rings."""
    smiles = "CCCCCCCCCCNC(=O)C[N+](C)(CCO)" \
             "CCCCCCCCCC[N+](C)(CCO)CC(=O)NCCCCCCCCCC"
    with pytest.raises(NovanaError):
        scaffold_from_smiles(smiles)


def test_exception_mixture_no_rings_shape():
    """Raises an error if no molecules with rings."""
    smiles = "OC(CC(=O)O)(CC(=O)O)C(=O)O.OC(CC(=O)O)(CC(=O)O)C(=O)O"
    with pytest.raises(NovanaError):
        shape_from_smiles(smiles)


def test_exception_molecule_no_rings_shape():
    """Raises an error if no molecules with rings."""
    smiles = "CCCCCCCCCCNC(=O)C[N+](C)(CCO)" \
             "CCCCCCCCCC[N+](C)(CCO)CC(=O)NCCCCCCCCCC"
    with pytest.raises(NovanaError):
        shape_from_smiles(smiles)
