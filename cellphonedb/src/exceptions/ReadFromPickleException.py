class ReadFromPickleException(Exception):
    def __init__(self, file='Can not read file'):
        super().__init__('Cannot parse pickle from file: {}'.format(file))
