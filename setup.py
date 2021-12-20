import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms_16S',
    version="0.0.2",
    description='OCMS_16S : Oxford Centre for Microbiome Studies pipelines for 16S amplicon sequencing analysis',
    author='Nicholas Ilott, Sandi Yen, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="microbiome, metagenomics, genomics",
    url="https://github.com/OxfordCMS/OCMS_16S",
    download_url="https://github.com/OxfordCMS/OCMS_16S/archive/refs/tags/v0.0.2.tar.gz",
    install_requires=["cgatcore"],
    packages=find_packages("./") + find_packages("./ocms16S/"),
    entry_points={
        'console_scripts': ['ocms_16s = ocms16S.ocms_16s:main']
    },
    include_package_data=True,
    python_requires='>=3.6.0'                                            
)

