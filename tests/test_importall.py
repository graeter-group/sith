from pathlib import Path
from importlib import import_module


def test_import_all():
    source_dir = Path(__file__).parent / '..' / 'src' / 'sith'
    
    
    for a in source_dir.rglob('*.py'):
        print('file ', a)
        if '__init__' not in str(a):
            module = 'sith.' + '.'.join(str(a).split('/')[9:]).replace('.py', '')
            
            print('module: ', module)
            try:
                import_module(module)
            except ImportError as e:
                assert False, f"Failed to import numpy: {e}"
