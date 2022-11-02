class EmptyResultException(Exception):
    def __init__(self, description: str = None, hint: str = None):
        super(EmptyResultException, self).__init__('Result is empty')
        self.description = description
        self.hint = hint
