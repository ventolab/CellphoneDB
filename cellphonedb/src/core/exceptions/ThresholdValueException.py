class ThresholdValueException(Exception):
    def __init__(self, threshold_value):
        super(ThresholdValueException, self).__init__(
            'Threshold value ({}) is not valid. Accepted range: 0<=threshold<=1'.format(threshold_value))
