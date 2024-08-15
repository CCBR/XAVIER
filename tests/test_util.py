import os
import warnings
from xavier.src.xavier.util import xavier_base


def test_xavier_base():
    assert str(xavier_base("a", "b", "c")).endswith("/a/b/c")
