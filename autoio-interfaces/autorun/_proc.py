""" Process stuff
"""

import os
import errno
import signal
from functools import wraps


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    """ check if process has died
    """

    def decorator(func):
        def _handle_timeout():
            print(error_message)
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator
