class ParseCountsException(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(ParseCountsException, self).__init__('Invalid Counts data')
        self.description = description
        self.hint = hint
