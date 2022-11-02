class ParseMetaException(Exception):
    def __init__(self):
        super(ParseMetaException, self).__init__('Error parsing Meta data')
