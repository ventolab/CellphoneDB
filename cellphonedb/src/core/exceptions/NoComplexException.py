class NoComplexException(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(NoComplexException, self).__init__('No Complexes data in the database')
        self.description = description
        self.hint = hint
