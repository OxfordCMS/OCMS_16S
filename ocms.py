'''
ocms.py - Oxford Centre for Microbiome Studies apps
=====================================================

'''

import os
import sys
import re
import glob
import imp
import cgatcore.iotools as iotools
import cgat


def main(argv=None):

    argv = sys.argv

    path = os.path.join(os.path.abspath(os.path.dirname(__file__)))

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print((globals()["__doc__"]))

    command = argv[1]
    #command = re.sub("-", "_", command)
    
    (file, pathname, description) = imp.find_module(command, [path])
    module = imp.load_module(command, file, pathname, description)

    # remove 'ocms' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
