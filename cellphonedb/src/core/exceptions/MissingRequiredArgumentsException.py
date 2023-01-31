class MissingRequiredArgumentsException(Exception):
    def __init__(self, description: str = None):
        super(MissingRequiredArgumentsException, self).__init__(description)
