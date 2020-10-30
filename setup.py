import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms 16S',
    version="0.0",
    description='ocms : Oxford Centre for Microbiome Studies apps',
    author='Nicholas Ilott, Sandi Yen, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="microbiome, metagenomics, genomics",
    packages=find_packages() + find_packages("./pipelines"),
    entry_points={
        'console_scripts': ['ocms = ocms:main']
    }
)
