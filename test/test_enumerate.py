from novana.api import enumerate_substructures_from_smiles
from utils import mol_list_to_smiles_list


def test_enumerate_kekule():
    smiles = "ClC1=CC2=C(C=C1)N(C(=O)N2)C"
    result = mol_list_to_smiles_list(
        enumerate_substructures_from_smiles(smiles))
    expected = ['Cn1c(=O)[nH]c2ccccc21', 'CN1CNc2cc(Cl)ccc21', 'O=c1[nH]c2ccc(Cl)cc2[nH]1']
    assert result == expected
