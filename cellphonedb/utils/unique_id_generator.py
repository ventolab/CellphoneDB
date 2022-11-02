import hashlib

import pandas as pd


def interaction(interaction: pd.Series, suffixes=('_x', '_y')) -> str:
    prefix = 'CPI'

    simple_complex_indicator_1 = _complex_indicator(interaction['is_complex{}'.format(suffixes[0])])
    simple_complex_indicator_2 = _complex_indicator(interaction['is_complex{}'.format(suffixes[1])])

    hash = _calculate_interaction_hash(interaction['name{}'.format(suffixes[0])],
                                       interaction['name{}'.format(suffixes[1])])

    control_digit = '0'

    return '{}-{}{}{}{}'.format(prefix, simple_complex_indicator_1, simple_complex_indicator_2, control_digit,
                                hash.upper())


def _complex_indicator(is_complex: bool) -> str:
    if is_complex:
        return 'C'
    else:
        return 'S'


def _calculate_interaction_hash(element_1: str, element_2: str, digits: int = 8) -> str:
    return hashlib.sha1((element_1 + element_2).encode('utf-8')).hexdigest()[:digits]
