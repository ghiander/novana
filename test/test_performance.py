import time
from novana.api import scaffold_from_smiles
from novana.api import shape_from_smiles


def time_function(func, *args, **kwargs):
    """Returns the time (ms) needed to run a function."""
    start = time.time()
    func(*args, **kwargs)
    end = time.time()
    return (end - start)*1000


def test_time_scaffold_2_rings():
    """Scaffold decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"
    t = time_function(scaffold_from_smiles, smiles)
    assert t < 2


def test_time_scaffold_4_rings():
    """Scaffold decomposition within a time cap."""
    smiles = "S2c1c(cc(cc1)C(=O)C)N(c3c2cccc3)CCCN4CCN(CC4)CCO"
    t = time_function(scaffold_from_smiles, smiles)
    assert t < 3


def test_time_scaffold_8_rings():
    """Scaffold decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"*4
    t = time_function(scaffold_from_smiles, smiles)
    assert t < 4


def test_time_scaffold_20_rings():
    """Scaffold decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"*10
    t = time_function(scaffold_from_smiles, smiles)
    assert t < 8


def test_time_shape_2_rings():
    """Shape decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"
    t = time_function(shape_from_smiles, smiles)
    assert t < 2


def test_time_shape_4_rings():
    """Shape decomposition within a time cap."""
    smiles = "S2c1c(cc(cc1)C(=O)C)N(c3c2cccc3)CCCN4CCN(CC4)CCO"
    t = time_function(shape_from_smiles, smiles)
    assert t < 3


def test_time_shape_8_rings():
    """Shape decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"*4
    t = time_function(shape_from_smiles, smiles)
    assert t < 4


def test_time_shape_20_rings():
    """Shape decomposition within a time cap."""
    smiles = "CC[n+]1ccccc1CNC2cccc2"*10
    t = time_function(shape_from_smiles, smiles)
    assert t < 8
