from novana.core import create_rwmol_from_smiles
from utils import mol_to_smiles


def test_largest_fragment_1():
    """The correct fragment is extracted from the mixture."""
    smiles = "N1CCNCC1.N2CCNCC2.N3CCNCC3.OC(CC(=O)O)" \
             "(CC(=O)O)C(=O)O.OC(CC(=O)O)(CC(=O)O)C(=O)O"
    result = mol_to_smiles(create_rwmol_from_smiles(smiles))
    expected = "C1CNCCN1"
    assert result == expected


def test_largest_fragment_2():
    """The correct fragment is extracted from the mixture."""
    smiles = "C=C1CN(S(=O)(=O)c2ccc(C)cc2)CCCN" \
             "(Cc2ccccc2)CCCN(S(=O)(=O)c2cccc(N(C)C)c2)C1.Cl"
    result = mol_to_smiles(create_rwmol_from_smiles(smiles))
    expected = "C=C1CN(S(=O)(=O)c2ccc(C)cc2)CCCN" \
               "(Cc2ccccc2)CCCN(S(=O)(=O)c2cccc(N(C)C)c2)C1"
    assert result == expected


def test_largest_fragment_3():
    """The correct fragment is extracted from the mixture."""
    smiles = "[s+]1c3c(nc2c1cc(cc2)N(C)C)ccc(c3)N(C)C" "." \
             "[N+](=O)([O-])c4c(ccc(c4)[N+](=O)[O-])" \
             "NN=C(c5ccc(cc5)O)c6ccc(cc6)O"
    result = mol_to_smiles(create_rwmol_from_smiles(smiles))
    expected = "O=[N+]([O-])c1ccc(NN=C(c2ccc(O)cc2)" \
               "c2ccc(O)cc2)c([N+](=O)[O-])c1"
    assert result == expected
