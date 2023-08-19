import os
import logging
from rdkit import RDLogger


def get_logger(name):
    """
    Creates a logger where its level is
    determined by the variable NOVANA_LOG_LEVEL.

    """
    log = logging.getLogger(name)
    if os.environ.get("NOVANA_LOG_LEVEL") == "DEBUG":
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return log


def disable_rdkit_logging():
    """Disables INFO and WARNING logs from RDKit."""
    RDLogger.DisableLog('rdApp.info')
    RDLogger.DisableLog('rdApp.warning')
