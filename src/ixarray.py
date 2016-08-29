#!/usr/bin/env python
import argparse
import sys
import IPython
import time
import xarray 

def initArgs():
    """Initialize argparse arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, 
                        help='Path')
    parser.add_argument("-f", "--file_name", type=str, 
                        help='File name')
    parser.add_argument("-i", "--instrument", type=str, 
                        help='Instrument')
    parser.add_argument("-e", "--exp", type=str, 
                        help='Experiment number')
    parser.add_argument("-r", "--run", type=int,  
                        help='Run number')
    parser.add_argument("--base", type=str, default='x', 
                        help='Base into which DataSource object is initiated.')
    return parser.parse_args()

def open_h5netcdf(file_name=None, path='', instrument=None, exp=None, run=None, **kwargs):
    if exp:
        if not instrument:
            instrument = exp[0:3]
    
    if not path:
        path = '/reg/d/psdm/{:}/{:}/scratch/nc/'.format(instrument, exp)

    if not file_name:
        file_name = '{:}run{:04}*.nc'.format(path, int(run))

    print file_name
    return xarray.open_mfdataset(file_name, engine='h5netcdf')

def banner(x, base='x', time0=None):
    print x
    print '='*80
    print 'See http://xarray.pydata.org for details on how to use data:'
    print '='*80

def main():
#    from pylab import *
    import pandas as pd
    time0 = time.time()
    args = initArgs()
    try:
        x = open_h5netcdf(**vars(args)) 
    except:
        import psxarray
        x = psxarray.get_xdat(**vars(args))

    setattr(sys.modules['__main__'], args.base, x)
    
    banner(x, base=args.base, time0=time0)
    #IPython.embed(banner1=banner(x, base=args.base, time0=time0))

if __name__ == "__main__":
    sys.exit(main())

