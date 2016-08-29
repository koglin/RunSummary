__all__ = []

import sys
sys.path = [xx for xx in sys.path if xx.find('koglin')>0] + sys.path
#sys.path = [xx for xx in sys.path if xx.find('xarrayenv')>0] + sys.path
#sys.path = [xx for xx in sys.path if xx.find('xarrayenv')>0 and xx.find('pandas')>0] + sys.path

from RunSummary import psxarray
from RunSummary import build_html

open_h5netcdf = psxarray.open_h5netcdf
Build_html = build_html.Build_html


__version__ = '00.00.01'

import logging

logger = logging.getLogger('PyDataSource')
logger.setLevel(logging.DEBUG)

#fh = logging.FileHandler('data_summary.log')
#fh.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
#fh.setFormatter(formatter)
ch.setFormatter(formatter)

#logger.addHandler(fh)
logger.addHandler(ch)


def set_logger_level(lvl):
    logger.setLevel( getattr(logging,lvl) )
#    fh.setLevel( getattr(logging,lvl) )
    ch.setLevel( getattr(logging,lvl) )
    return

def logger_flush():
#    fh.flush()
    return

def initArgs():
    """Initialize argparse arguments.
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--exp", type=str, 
                        help='Experiment number')
    parser.add_argument("-r", "--run", type=int,  
                        help='Run number')
    parser.add_argument("--nchunks", type=int,  
                        help='total number of chunks')
    parser.add_argument("--ichunk", type=int,  
                        help='chunk index')
    parser.add_argument("-i", "--instrument", type=str, 
                        help='Instrument')
    parser.add_argument("-s", "--station", type=int, 
                        help='Station')
    parser.add_argument("-n", "--nevents", type=int, 
                        help='Number of events to analyze')
    parser.add_argument("--make_summary", action="store_true", default=False,
                        help='Make summary for array data.')
    parser.add_argument("--show_errors", action="store_true", default=False,
                        help='Show Errors in cases that might not be explicit ' \
                             'due to try/except statements')
    return parser.parse_args()



