class NoInteractionsFound(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(NoInteractionsFound, self).__init__('No CellphoneDB interactions found in this input.')
        self.description = description
        self.hint = hint
