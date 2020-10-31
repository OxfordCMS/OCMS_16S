import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms_16S',
    version="0.0.1",
    description='OCMS_16S : Oxford Centre for Microbiome Studies pipelines for 16S amplicon sequencing analysis',
    author='Nicholas Ilott, Sandi Yen, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="microbiome, metagenomics, genomics",
    url="https://github.com/OxfordCMS/OCMS_16S",
    packages=find_packages("./") + find_packages("./pipelines"),
    entry_points={
        'console_scripts': ['ocmsflow = ocms_16S.ocmsflow:main']
    },
    include_package_data=True,
    python_requires='>=3.6.0'                                            
)

