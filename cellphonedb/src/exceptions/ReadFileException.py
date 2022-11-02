class ReadFileException(Exception):
    def __init__(self, file='Can not read file'):
        super(ReadFileException, self).__init__('Can not read {}'.format(file))
