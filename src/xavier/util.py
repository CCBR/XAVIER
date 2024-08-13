#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import sys
import hashlib
import subprocess
import json
import glob
import os
import warnings


def xavier_base(*paths):
    """Get the absolute path to a file in the repository
    @return abs_path <str>
    """
    basedir = pathlib.Path(__file__).absolute().parent.parent.parent
    return basedir.joinpath(*paths)


def get_version():
    """Get the current version
    @return version <str>
    """
    with open(xavier_base("VERSION"), "r") as vfile:
        version = f"v{vfile.read().strip()}"
    return version
