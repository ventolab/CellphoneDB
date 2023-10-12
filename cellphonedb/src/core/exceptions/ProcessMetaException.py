class ProcessMetaException(Exception):
    def __init__(self):
        super(ProcessMetaException, self).__init__(
            'Error processing Meta data - please consult ' +
            'https://cellphonedb.readthedocs.io/en/latest/RESULTS-DOCUMENTATION.html#meta-file')
