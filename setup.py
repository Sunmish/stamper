#! /usr/bin/env python

import setuptools
import sys

with open("README.md", "r") as f:
    long_description = f.read()

if sys.version_info[0] == 3:
    _version = ">=3.0.0"
elif sys.version_info[0] == 2:
    _version = "<3.0.0"
else:
    _version = ""

reqs = [
    "astropy"+_version,
    "numpy",
    "mechanize",
]


setuptools.setup(
    name="stamper",
    version="0.0.0",
    author="Stefan W Duchesne",
    author_email="stefanduchesne@gmail.com",
    description="A small package to download postage stamps.",
    long_description=long_description,
    url="https://github.com/Sunmish/stamper",
    install_requires=reqs,
    packages=["stamper"],
)
