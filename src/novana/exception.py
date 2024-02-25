class NovanaError(Exception):
    pass


class NonCyclicMoleculeError(Exception):
    def __init__(self, message="The input does not contain "
                               "any cyclic molecule."):
        super().__init__(message)
