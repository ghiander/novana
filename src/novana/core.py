from novana.chem import MoleculeBuilder
from novana.process import TerminalAtomRemover
from novana.process import SingleBondedCarborConverter
from novana.valence import HydrogenAdjuster
from novana.utils import get_logger
from novana.utils import disable_rdkit_logging


# Configure logging
log = get_logger(__name__)
disable_rdkit_logging()


def create_rwmol_from_smiles(smiles, flatten_mol=True):
    return MoleculeBuilder.create_rwmol_from_smiles(smiles, flatten_mol)


def adjust_hydrogens(rwmol):
    """Adjusts hydrogens in the molecule using heuristics."""
    return HydrogenAdjuster(rwmol).process_and_return()


def recursively_remove_single_bonded_atoms(rwmol):
    """Removes recursively all single-bonded atoms and their bonds."""
    return TerminalAtomRemover(rwmol).process_and_return()


def convert_into_single_bonded_carbons(rwmol):
    """Converts all atoms into single-bonded carbons."""
    return SingleBondedCarborConverter(rwmol).process_and_return()
