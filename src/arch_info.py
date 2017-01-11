from pylab import *
import pandas as pd

def load_archive_pvs(ds, pvs=None, quiet=True):
    """load archive data into Dataset
    """
    import xarray as xr
    meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'}
    arch_data = get_pvs(ds, pvs=pvs, quiet=quiet)
    data_arrays = []
    data_points = {}
    for alias, dat in arch_data.items():
        attrs = {a: dat['meta'].get(val) for a,val in meta_attrs.items() if val in dat['meta']}
        vals = [item['val'] for item in dat['data']]
        if isinstance(vals[0],str):
            if not quiet:
                print alias, 'string'
            vals = np.array(vals, dtype=str)
        else:
            times = [np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
            times = np.array(times, dtype=times[0].dtype)
            if len(vals) > 1:
                data_arrays.append( xr.DataArray(vals, coords=[times], dims=['time'], name=alias, attrs=attrs) )
            else:
                data_points[alias] = {'value': vals[0], 'time': times[0], 'pv': attrs['pv'], 'units': attrs.get('units','')}

    pv_dset = xr.merge(data_arrays)
    pv_dset = pv_dset[sorted(pv_dset.data_vars)]
    #pv_sum = to_summary(data_arrays).to_dataframe().T
    #units = {attr: data_arrays[attr].attrs.get('units','') for attr in data_arrays}

    return pv_dset, pd.DataFrame(data_points).T #, pv_sum

#def get_pv_changes(
#[attr for attr in attrs if a[attr].dropna('time').values.std() < 10**(-int(a[attr].attrs.get('PREC',0)))]

def get_pvs(ds=None, start_time=None, end_time=None, pvs=None, quiet=True):
    """Get pvs from epics arch
    """
    from RunSummary.epicsarchive import EpicsArchive
    arch = EpicsArchive()
    arch_data = {}
    try:
        if not start_time:
            start_time = min(ds.scanData.start_times)
        if not end_time:
            end_time = max(ds.scanData.end_times)
    except:
        print '{:} Not a valid PyDataSource.DataSource'.format(ds)
        print 'Must supply both start_time and end_time or ds'

    tstart = sec_to_array(start_time)
    tend = sec_to_array(end_time)
    #tstart = np.datetime64(int(start_time), 's').tolist()
    #tend = np.datetime64(int(end_time), 's').tolist()

    # get pvs in epicsarchiver that are in data_source (i.e., in daq epicsArch.txt file)
    if not pvs:
        if not quiet:
            print 'Getting PVs from', ds
        pvs = _get_pv_dict(ds, quiet=quiet)
    
    # add set VAL for motor pvs.
    for pv, alias in pvs.items():
        if pv.endswith('.RBV'):
            pv_val = pv.rstrip('RBV')+'VAL'
            if arch.search_pvs(pv_val, do_print=False) != []: 
                pvs[pv_val] = alias+'_set'

    for pv, alias in pvs.items():
        if not quiet:
            print pv, alias
        try:
            arch_data[alias] = arch._get_json(pv, tstart, tend, False)[0] 
        except:
            if not quiet:
                print pv, alias, 'No Data'

    return arch_data

def _get_pv_dict(ds, quiet=True):
    """Get a dict of pvs from a PyDataSource.DataSource object
    """
    from RunSummary.epicsarchive import EpicsArchive
    arch = EpicsArchive()
    pvs = {pv:ds.epicsData.alias(pv) for pv in ds.epicsData.pvNames() if arch.search_pvs(pv, do_print=False) != []} 
    if not quiet:
        print '{:} contains {:} pvs'.format(ds, len(pvs))

    return pvs

def sec_to_array(sec):
    """Convert datetime64 to time list.
    """
    dt = np.datetime64(int(sec), 's').tolist()
    return [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]


