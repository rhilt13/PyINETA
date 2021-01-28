"""
    pyINETA for Python3
    Author : Rahil Taujale
"""
import sys

__release__ = '201028'
__version__ = '2.0'

sys.stderr.write("""PyINETA [ Release %s ]
import sys 
sys.dont_write_bytecode = True
import numpy as np
np.set_printoptions(suppress=True)
Dependencies (update as i go along?):
sudo python3 -m pip install numpy pandas scipy bokeh panel pillow
""" % __version__)





# from .read import *

# from .alignment.constructor import AlignmentArray as AlignmentArray
# from .read.xma import xmaReader as xma


# from .fig.comparitor import CompareLogo


# # no integration
# from .seqnames import SeqNames
# from .taxdump import taxdump
# from .pdbmeta import pdbmeta