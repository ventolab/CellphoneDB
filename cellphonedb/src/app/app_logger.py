import logging

import sys

app_logger = logging.getLogger(__name__)


def setLevel(level: str = 'INFO'):
    app_logger.setLevel(getattr(logging, level))


if not app_logger.handlers:
    app_logger.propagate = False
    ch = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('[ ][APP][%(asctime)s][%(levelname)s] %(message)s', "%d/%m/%y-%H:%M:%S")
    ch.setFormatter(formatter)

    setLevel()
    app_logger.addHandler(ch)
