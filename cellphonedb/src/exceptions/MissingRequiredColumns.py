class MissingRequiredColumns(Exception):

    def __init__(self, columns: list):
        super().__init__('Missing columns: `{}`'.format(','.join(columns)))