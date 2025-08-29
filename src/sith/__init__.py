from pathlib import Path


class variables4tests:
    resources = Path(__file__).parent / '..' / '..' / 'tests' / 'resources'
    G = resources / 'G'