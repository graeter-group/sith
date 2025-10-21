from .SITH import SITH
from pathlib import Path


SITH.__module__ = __name__

class variables4tests:
    resources = Path(__file__).parent / '..' / '..' / 'tests' / 'resources'
    G = resources / 'G'
