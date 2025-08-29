from pathlib import Path
from importlib import import_module
import importlib.util as iu

source_dir = Path(__file__).parent / '..'

for a in source_dir.rglob('*.py'):
    if '__init__' not in str(a):
        module = 'sith.' + '.'.join(str(a).split('/')[9:]).replace('.py', '')
        if '.tests.' not in module:
            print(module)
            try:
                import_module(module)
            except ImportError as e:
                assert False, f"Failed to import numpy: {e}"
