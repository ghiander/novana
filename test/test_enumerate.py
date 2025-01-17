from novana.api import substructure_scaffolds_from_smiles
from utils import mol_list_to_smiles_list


def test_enumerate_kekule():
    smiles = "ClC1=CC2=C(C=C1)N(C(=O)N2)C"
    result = mol_list_to_smiles_list(
        substructure_scaffolds_from_smiles(smiles))
    expected = ['Cn1c(=O)[nH]c2ccccc21', 'CN1CNc2cc(Cl)ccc21', 'O=c1[nH]c2ccc(Cl)cc2[nH]1']
    assert result == expected


def test_enumerate_remove_linker_product():
    smiles = "C(C)(C)COC1CC[C@H]2CN(C[C@H]2C1)C(=O)C[C@H](CC3=C(F)C=CC=C3Br)C(=O)O"
    result = mol_list_to_smiles_list(
        substructure_scaffolds_from_smiles(smiles))
    expected = ['O=C(O)C(CC(=O)N1CC2CCCCC2C1)Cc1c(F)cccc1Br',
                'CC(C)COC1CCC2CN(C(=O)CC(Cc3ccccc3F)C(=O)O)CC2C1',
                'CC(C)COC1CCC2CN(C(=O)CC(Cc3ccccc3Br)C(=O)O)CC2C1']
    assert result == expected
