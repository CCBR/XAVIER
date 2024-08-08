import os
import warnings
from xavier.src.xavier.util import (
    xavier_base
)

def test_xavier_base():
    test_base = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    assert all([xavier_base() == test_base,
                xavier_base("a","b","c") == os.path.join(test_base, "a", "b", "c")
                ])
