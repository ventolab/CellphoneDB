class AllCountsFilteredException(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(AllCountsFilteredException, self).__init__('All counts filtered')
        self.description = description
        self.hint = hint
