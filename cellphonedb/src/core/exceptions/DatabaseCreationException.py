class DatabaseCreationException(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(DatabaseCreationException, self).__init__(
            "Failed to create CellphoneDB database due to failed data sanity tests")
        self.description = description
        self.hint = hint
