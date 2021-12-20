## This script was lifted from https://github.com/cgat-developers/cgat-apps
## in order to remove the cgat-apps dependency from OCMS_16S and without
## having to re-write the tool.

## Copyright (c) 2021 cgat-developers

'''combine_tables.py - join tables
===============================

:Tags: Python

Purpose
-------

This script reads several tab-separated tables and joins them into a
single one.

.. todo::

   * Rename to tables2table.py
   * Use pandas dataframes for fast IO and merging/joining

Usage
-----

The option ``--header-names`` sets the column titles explicitely. Add
``--skip-titles`` if you want to avoid echoing the original title in
the input files.


Example::

   python combine_tables.py --help

Type::

   python combine_tables.py --help

for command line help.

Command line options
--------------------

'''

import sys
from cgatcore.tables import main


if __name__ == '__main__':
    sys.exit(main(sys.argv))
