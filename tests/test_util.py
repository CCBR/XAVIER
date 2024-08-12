import os
import warnings
from xavier.src.xavier.util import (
    xavier_base
)

def test_xavier_base():
    test_base = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    xavier_base()
    assert xavier_base("a","b","c").endswith('/a/b/c')
