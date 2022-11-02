class NoInteractionsFound(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(NoInteractionsFound, self).__init__('No CellPhoneDB interacions found in this input.')
        self.description = description
        self.hint = hint
