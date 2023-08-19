from novana.api import scaffold_from_smiles
from novana.api import shape_from_smiles
from utils import mol_to_inchikey
from utils import mol_to_smiles


def test_scaffold_1():
    """Scaffold decomposition."""
    smiles = "COC(=O)N1CC2(C1)CS(=O)(=Nc1cc(C)c3c" \
             "(Nc4ccc(F)cc4O[C@H](C)C(=O)NCC(F)(F)F)ncnc3c1)C2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(Nc2ncnc3cc(N=S4CC5(CNC5)C4)ccc23)cc1"
    assert result == expected


def test_scaffold_2():
    """Scaffold decomposition."""
    smiles = "P(=O)(OCc1ccc(cc1)OC(=O)c2ccccc2)" \
             "(OCc3ccc(cc3)OC(=O)c4ccccc4)CCCC(=O)N(O)C"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(COc2ccc(COPOCc3ccc(OCc4ccccc4)cc3)cc2)cc1"
    assert result == expected


def test_scaffold_3():
    """Scaffold decomposition."""
    smiles = "C=C1CN(S(=O)(=O)c2ccc(C)cc2)CCCN" \
             "(Cc2ccccc2)CCCN(S(=O)(=O)c2cccc(N(C)C)c2)C1"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(CN2CCCN(Sc3ccccc3)CCCN(Sc3ccccc3)CCC2)cc1"
    assert result == expected


def test_scaffold_4():
    """Scaffold decomposition."""
    smiles = "[O-]S(=O)(=O)C(O)[C@@H](NC(=O)[C@H](NC(=O)" \
             "OC1(CCN(CC1)C(=O)OC(C)(C)C)CC)CC(C)C)C[C@@H]2CCNC2=O"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C(CNCOC1CCNCC1)NCCC1CCNC1"
    assert result == expected


def test_scaffold_5():
    """Scaffold decomposition."""
    smiles = "C[C@@H]1CC[C@@]23CCC(=O)[C@H]2[C@@]1(C)" \
             "[C@@H](C[C@@](C)(C=C)[C@@H](O)[C@@H]3C)OC=O"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1CCC2CCCC3(CC1)CCCC23"
    assert result == expected


def test_scaffold_6():
    """Scaffold decomposition."""
    smiles = r"CC1=CCC(O)/C=C\C(C)C(O)C(C)C=C(C)C(=O)" \
             r"c2c(O)c(C)cc3c2C(=O)C(N)=C(NC(=O)" \
             r"/C=C\C=C/C=C\C(C)C(O)CC1=O)C3=O"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = r"C1=CCCCC/C=C\C=C/C=C\CNC2=CCc3c(cccc3C2)CC=CCCC/C=C\CC1"
    assert result == expected


def test_scaffold_7():
    """Scaffold decomposition."""
    smiles = "CCN(CC)[c+]1sc2c(s1)C(c1ccc(Cl)cc1)Oc1c(I)cc(I)cc1-2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(C2Oc3ccccc3C3=C2SCS3)cc1"
    assert result == expected


def test_scaffold_8():
    """Scaffold decomposition."""
    smiles = "COC(=O)N1CC2(C1)CS(=O)(=Nc1cc(C)c3c(Nc4ccc(F)" \
             "cc4O[C@H](C)C(=O)NCC(F)(F)F)ncnc3c1)C2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(Nc2ncnc3cc(N=S4CC5(CNC5)C4)ccc23)cc1"
    assert result == expected


def test_scaffold_9():
    """Scaffold decomposition."""
    smiles = "P(=O)(OCc1ccc(cc1)OC(=O)c2ccccc2)(OCc3ccc" \
             "(cc3)OC(=O)c4ccccc4)CCCC(=O)N(O)C"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(COc2ccc(COPOCc3ccc(OCc4ccccc4)cc3)cc2)cc1"
    assert result == expected


def test_scaffold_mixture_1():
    """Scaffold decomposition."""
    smiles = "S2c1c(cc(cc1)C(=O)C)N(c3c2cccc3)CCCN4CCN(CC4)CCO" "." \
             r"OC(=O)/C=C\C(=O)O" "." r"OC(=O)/C=C\C(=O)O"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc2c(c1)Sc1ccccc1N2CCCN1CCNCC1"
    assert result == expected


def test_scaffold_mixture_2():
    """Scaffold decomposition."""
    smiles = "C=C1CN(S(=O)(=O)c2ccc(C)cc2)CCCN(Cc2ccccc2)" \
             "CCCN(S(=O)(=O)c2cccc(N(C)C)c2)C1" "." "Cl"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "c1ccc(CN2CCCN(Sc3ccccc3)CCCN(Sc3ccccc3)CCC2)cc1"
    assert result == expected


def test_scaffold_mixture_3():
    """Scaffold decomposition."""
    smiles = "[Na+].[O-]S(=O)(=O)C(O)[C@@H](NC(=O)[C@H](NC(=O)" \
             "OC1(CCN(CC1)C(=O)OC(C)(C)C)CC)CC(C)C)C[C@@H]2CCNC2=O"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C(CNCOC1CCNCC1)NCCC1CCNC1"
    assert result == expected


def test_scaffold_mixture_4():
    """Scaffold decomposition."""
    smiles = "[I-].CC[n+]1ccccc1CNC2cccc2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1=CC(NCc2ccccn2)C=C1"
    assert result == expected


def test_scaffold_charged_1():
    """Scaffold decomposition."""
    smiles = "[BH2-]1n2cccc2C(c2ccccc2)=C2C=CC=[N+]12"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "B1n2cccc2C(c2ccccc2)=C2C=CC=[N+]12"
    assert result == expected


def test_scaffold_charged_2():
    """Scaffold decomposition."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1=CC(NCc2ccccn2)C=C1"
    assert result == expected


def test_scaffold_charged_3():
    """Scaffold decomposition."""
    smiles = "CCN(CC)[c+]1sccs1"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1=CSCS1"
    assert result == expected


def test_scaffold_charged_4():
    """Scaffold decomposition."""
    smiles = "C[C@H]1[C@H]2[C@H](C[C@H]3[C@@H]4CC=C5C[C@H]" \
             "(n6cc[n+](Cc7ccccc7Br)c6)CC[C@]5(C)" \
             "[C@H]4CC[C@@]32C)O[C@]12CC[C@@H](C)CO2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1=C2CC(n3cc[n+](Cc4ccccc4)c3)" \
               "CCC2C2CCC3C4CC5(CCCCO5)OC4CC3C2C1"
    assert result == expected


def test_scaffold_charged_5():
    """Scaffold decomposition."""
    smiles = "[BH2-]1n2cccc2C(c2ccccc2)=C2C=CC=[N+]12"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "B1n2cccc2C(c2ccccc2)=C2C=CC=[N+]12"
    assert result == expected


def test_scaffold_charged_6():
    """Scaffold decomposition."""
    smiles = "C[C@H]1[C@H]2[C@H](C[C@H]3[C@@H]4CC=C5C[C@H]" \
             "(n6cc[n+](Cc7ccccc7Br)c6)CC[C@]5(C)[C@H]4CC" \
             "[C@@]32C)O[C@]12CC[C@@H](C)CO2"
    result = mol_to_smiles(scaffold_from_smiles(smiles))
    expected = "C1=C2CC(n3cc[n+](Cc4ccccc4)c3)" \
               "CCC2C2CCC3C4CC5(CCCCO5)OC4CC3C2C1"
    assert result == expected


def test_shape_1():
    """Shape decomposition."""
    smiles = "CCN(CC)[c+]1sc2c(s1)C(c1ccc(Cl)cc1)Oc1c(I)cc(I)cc1-2"
    result = mol_to_smiles(shape_from_smiles(smiles))
    expected = "C1CC2C3CCCCC3CC(C3CCCCC3)C2C1"
    assert result == expected


def test_shape_2():
    """Shape decomposition."""
    smiles = "C[C@@H]1CC[C@@]23CCC(=O)[C@H]2[C@@]1(C)" \
             "[C@@H](C[C@@](C)(C=C)[C@@H](O)[C@@H]3C)OC=O"
    result = mol_to_smiles(shape_from_smiles(smiles))
    expected = "C1CCCC2CCCC3(C1)CCCC23"
    assert result == expected


def test_shape_3():
    """Shape decomposition."""
    smiles = r"CC1=CCC(O)/C=C\C(C)C(O)C(C)C=C(C)C(=O)" \
             r"c2c(O)c(C)cc3c2C(=O)C(N)=C(NC(=O)" \
             r"/C=C\C=C/C=C\C(C)C(O)CC1=O)C3=O"
    result = mol_to_smiles(shape_from_smiles(smiles))
    expected = "C1CCCCCCCCCCCCC2CCC3C(CCCCCCCCCCC1)CCCC3C2"
    assert result == expected


def test_identical_scaffold_shape_1():
    """Scaffold and shape are identical."""
    smiles = "C[C@@H]1CC[C@@]23CCC(=O)[C@H]2[C@@]1(C)" \
             "[C@@H](C[C@@](C)(C=C)[C@@H](O)[C@@H]3C)OC=O"
    scaffold = mol_to_inchikey(scaffold_from_smiles(smiles))
    shape = mol_to_inchikey(shape_from_smiles(smiles))
    assert scaffold == shape


def test_different_scaffold_tautomers_1():
    """Tautomers yield different scaffolds."""
    smiles = r"CC1=CCC(O)/C=C\C(C)C(O)C(C)C=C(C)C(=O)" \
             r"c2c(O)c(C)cc3c2C(=O)C(N)=C(NC(=O)" \
             r"/C=C\C=C/C=C\C(C)C(O)CC1=O)C3=O"
    result_1 = mol_to_inchikey(scaffold_from_smiles(smiles))
    smiles = r"CC1\C=C/C=C\C=C/C(=NC2=C(N)C(=O)c3c(cc(C)c" \
             r"(O)c3C(=O)C(=CC(C)C(O)C(C)\C=C/C(O)CC=C(C)" \
             r"C(=O)CC1O)C)C2=O)O"
    result_2 = mol_to_inchikey(scaffold_from_smiles(smiles))
    assert result_1 != result_2


def test_identical_shape_tautomers_1():
    """Tautomers yield same shape."""
    smiles = r"CC1=CCC(O)/C=C\C(C)C(O)C(C)C=C(C)C(=O)" \
             r"c2c(O)c(C)cc3c2C(=O)C(N)=C(NC(=O)" \
             r"/C=C\C=C/C=C\C(C)C(O)CC1=O)C3=O"
    result_1 = mol_to_inchikey(shape_from_smiles(smiles))
    smiles = r"CC1\C=C/C=C\C=C/C(=NC2=C(N)C(=O)c3c(cc(C)c" \
             r"(O)c3C(=O)C(=CC(C)C(O)C(C)\C=C/C(O)CC=C(C)" \
             r"C(=O)CC1O)C)C2=O)O"
    result_2 = mol_to_inchikey(shape_from_smiles(smiles))
    assert result_1 == result_2
