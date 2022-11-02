class NotADataFrameException(Exception):
    def __init__(self, file='Can not read file'):
        super().__init__('Given file doesn\'t contain a pandas DataFrame: {}'.format(file))
