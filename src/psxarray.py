# standard python modules
import os
import operator
import time
import cPickle as pickle

# scipy and xarray
from pylab import *
import xarray

#PyDatasource modules
import PyDataSource
from welford import Welford

#psana methods
from PSCalib.GeometryAccess import GeometryAccess, img_from_pixel_arrays

#ds = PyDataSource.DataSource('exp=cxij4915:run=49:smd')

def load_file(file_name=None, run=None, exp=None, path=None):
    if not path:
        if exp:
            instrument = exp[0:3]
            path = '/reg/d/psdm/{:}/{:}/scratch/nc/'.format(instrument,exp)
        else:
            path = ''
                
    if not file_name:
        if not run:
            print 'file_name or run number must be supplied.'
            return 
        else:
            file_name = 'run{:04}.pkl'.format(run)
   
    print path+file_name
    with open(path+file_name, 'rb') as pickle_file:
        data = pickle.load(pickle_file)

    return pickle.loads(data)

def write_file(data, file_name=None, run=None, exp=None, path=None):
    if not path:
        if not exp:
            exp = data.experiment

        if exp:
            instrument = exp[0:3]

        path = '/reg/d/psdm/{:}/{:}/scratch/'.format(instrument,exp)
                
    if not file_name:
        if not run:
            run = int(data.run[0])

        file_name = 'run{:04}.pkl'.format(run)
    
    with open(path+file_name, 'wb') as pickle_file:
        pickle.dump(pickle.dumps(data, protocol=-1), pickle_file, protocol=-1)

    return

def open_h5netcdf(file_name=None, path='', exp=None, run=None, chunk=False):
    if exp:
        instrument = exp[0:3]
    
    if not path:
        path = '/reg/d/psdm/{:}/{:}/scratch/nc/'.format(instrument, exp)

    if chunk and run:
        file_base = '{:}/run{:04}*.nc'.format(path, int(run))
        return xarray.open_mfdataset(file_base, engine='h5netcdf')
        
    else:
        if not file_name:
            file_name = '{:}/run{:04}.nc'.format(path, int(run))

        return xarray.open_dataset(file_name, engine='h5netcdf')

def to_h5netcdf(xdat=None, ds=None, file_name=None, path=None, **kwargs):
    if xdat:
        if xdat.__class__.__name__ == 'DataSource':
            # 1st arg is actually PyDataSource.DataSource
            ds = xdat
            xdat = None

    if not xdat:
        if not ds:
            ds = PyDataSource.DataSource(**kwargs)

        xdat = get_xdat(ds, **kwargs)
    
    if not path:
        path = '/reg/d/psdm/{:}/{:}/scratch/nc/'.format(xdat.instrument,xdat.experiment)

    if not os.path.isdir(path):
        os.mkdir(path)

    if not file_name:
        if 'ichunk' in xdat.attrs:
            file_name = '{:}/run{:04}_c{:03}.nc'.format(path, int(xdat.run[0]), xdat.attrs['ichunk'])

        else:
            file_name = '{:}/run{:04}.nc'.format(path, int(xdat.run[0]))
    
    xdat.to_netcdf(file_name, engine='h5netcdf')

def add_index(x, attr, name=None, nbins=8, bins=None, percentiles=None):

    if not bins:
        if not percentiles:
            percentiles = (arange(nbins+1))/float(nbins)*100.

        bins = np.percentile(x[attr].to_pandas().dropna(), percentiles)

    if not name:
        name = attr+'_index'

    x[name] = (['time'], np.digitize(x[attr].values, bins))


def get_xdat(ds=None, nevents=None, max_size=10001, 
        xbase=None, 
        publish=False,
        make_summary=False,
        ichunk=None,
        nchunks=24,
        #code_flags={'XrayOff': [162], 'XrayOn': [-162], 'LaserOn': [183, -162], 'LaserOff': [184, -162]},
        code_flags={'XrayOff': [162], 'XrayOn': [-162]},
        eventCodes=None, config=None, make_images=True, **kwargs):
    """Get xarray data.
       ichunk: chunk index (skip ahead nevents*ichunk)
    """
    if not ds:
        ds = PyDataSource.DataSource(**kwargs)

    if not nevents:
        if ichunk is not None:
            nevents = int(np.ceil(ds.nevents/float(nchunks)))
            print 'Do {:} of {:} events for {:} chunk'.format(nevents, ds.nevents, ichunk)
        else:
            nevents = ds.nevents
   
   

    adat = {}
    sumdat = {}
    ds.reload()
    evt = ds.events.next(publish=publish, init=publish)
    dtime = evt.EventId
    if not eventCodes:
        eventCodes = sorted(ds.configData._eventcodes.keys())
    
    if hasattr(ds.configData, 'ScanData') and ds.configData.ScanData:
        nsteps = ds.configData.ScanData.nsteps
    else:
        nsteps = 1
    
    neventCodes = len(eventCodes)
    det_funcs = {}

    for srcstr, src_info in ds.configData._sources.items():
        det = src_info['alias']
        nmaxevents = 100
        ievt = 0
        try:
            while det not in evt._attrs and ievt < nmaxevents:
                evt.next(publish=publish, init=publish)
                ievt += 1
            
            detector = getattr(evt,det)
            if hasattr(detector, '_update_xarray_info'):
                print 'updating', srcstr, det
                detector._update_xarray_info()

            # Note that the a and b objects here link the det_funcs dictionary 
            # to the adat dictionary of xarray.Dataset objects.
            # Thus, updating the det_funcs dictionary updates the xarray.Dataset objects.
            adat[det] = {}
            sumdat[det] = {}
            det_funcs[det] = {}
            xarray_dims = detector._xarray_info.get('dims')
            if xarray_dims is not None: 
                for attr,item in sorted(xarray_dims.items(), key=operator.itemgetter(0)):
                    # Only save data with less than max_size total elements
                    alias = det+'_'+attr
                    if len(item) == 3:
                        attr_info = item[2]
                    else:
                        attr_info = {}

                    det_funcs[det][attr] = {'alias': alias, 'det': det, 'attr': attr, 'attr_info': attr_info}
                    if np.product(item[1]) <= max_size:
                        a = [det+'_'+name for name in item[0]]
                        a.insert(0, 'time')
                        try:
                            b = list(item[1])
                        except:
                            b = [item[1]]
                        b.insert(0, nevents)

                        adat[det][alias] = (a, tuple(b))
                        det_funcs[det][attr]['event'] = {'dims': a, 'shape': b}

                    # Make avg and std objects for each eventCode and step for multi-dimension data
                    if make_summary and len(item[1]) > 1:
                        a = [det+'_'+name for name in item[0]]
                        a.insert(0, 'codes')
                        a.insert(0, 'steps')
                        #a.insert(0, 'codes')
                        #a.insert(0, 'steps')
                        b = list(item[1])
                        b.insert(0, neventCodes)
                        b.insert(0, nsteps)
                        summary_funcs = []
                        for istep in range(nsteps):
                            summary_funcs.append([Welford() for ec in eventCodes])
                        
                        summary_names = ['ec{:}'.format(ec) for ec in eventCodes]
                        det_funcs[det][attr]['welford'] = {'dims': a, 'shape': b, 
                                'names': summary_names,
                                'funcs': summary_funcs,
                                'attrs': ['mean', 'std', 'min', 'max']}
                        for fattr in det_funcs[det][attr]['welford']['attrs']:
                            sumdat[det][alias+'_'+fattr] = (a, tuple(b))

        except:
            print 'ERROR loading', srcstr, det

    #return sumdat, adat, det_funcs

    axdat = {}
    atimes = {}
    btimes = []
    if not xbase:
        xbase = xarray.Dataset()
    
    # Experiment Attributes
    xbase.attrs['data_source'] = str(ds.data_source)
    xbase.attrs['run'] = ds.data_source.run
    for attr in ['instrument', 'experiment', 'expNum', 'calibDir']:
        xbase.attrs[attr] = getattr(ds, attr)

    xbase.coords['time'] = np.zeros(nevents, dtype=dtime.datetime64.dtype)
    ttypes = {'sec': 'int32', 
              'nsec': 'int32', 
              'fiducials': 'int32', 
              'ticks': 'int32', 
              'run': 'int32'}
    # explicitly order EventId coords in desired order 
    #print nevents, ttypes
    #return xbase
    for attr in ['sec', 'nsec', 'fiducials', 'ticks', 'run']:
        dtyp = ttypes[attr]
        print attr, dtyp
        xbase.coords[attr] = (['time'], np.zeros(nevents,dtype=dtyp))

    xbase.coords['step'] = (['time'], np.empty(nevents,dtype=int))
    
    # Event Codes -- unfortunately cannot use bool in netcdf4 so use byte
    for code in eventCodes:
        xbase.coords['ec{:}'.format(code)] = ('time', np.zeros(nevents, dtype=byte))

    for attr, ec in code_flags.items():
        xbase.coords[attr] = ('time', np.zeros(nevents, dtype=byte))
        xbase.coords[attr].attrs['doc'] = 'Event code flag: True if all positive and no negative "codes" are in eventCodes'
        xbase.coords[attr].attrs['codes'] = ec

    xbase.attrs['event_flags'] = code_flags.keys()

    xbase.coords['steps'] = range(nsteps)
    xbase.coords['codes'] = eventCodes
 
    # Scan Attributes -- cannot put None or dicts as attrs in netcdf4
    # e.g., pvAliases is a dict
    if hasattr(ds.configData, 'ScanData') and ds.configData.ScanData:
        if ds.configData.ScanData.nsteps == 1:
            attrs = ['nsteps']
        else:
            attrs = ['nsteps', 'pvControls', 'pvMonitors', 'pvLabels']
            #xbase.coords['pvControls'] = ds.configData.ScanData.pvControls
            for attr, vals in ds.configData.ScanData.control_values.items():
                alias = ds.configData.ScanData.pvAliases[attr]
                xbase.coords[alias+'_steps'] = (['steps'], vals) 
                xbase.coords[alias] = ('time', np.zeros(nevents, dtype=byte))

        for attr in attrs:
            val = getattr(ds.configData.ScanData, attr)
            if val:
                xbase.attrs[attr] = val 

   
    for srcstr, item in sorted(ds.configData._sources.items(), key=operator.itemgetter(0)):
        det = item['alias']
        if det in adat:
            atimes[det] = []
            axdat[det] = xarray.Dataset()
            for attr, item in sorted(adat[det].items(), key=operator.itemgetter(0)):
                axdat[det][attr] = (item[0], np.zeros(item[1]))
            
            axdat[det].coords['steps'] = range(nsteps)
            axdat[det].coords['codes'] = eventCodes

    for srcstr, srcitem in sorted(ds.configData._sources.items(), key=operator.itemgetter(0)):
        src_info = ds.configData._sources.get(srcstr).copy()
        src_info['src'] = str(src_info.get('src'))
        det = srcitem.get('alias')
        readoutCode = srcitem.get('eventCode')
        detector = getattr(evt, det)
        det_func = det_funcs.get(det)
        if det_func and det in adat:
            while det not in evt._attrs:
                evt.next(publish=publish, init=publish)

            detector._update_xarray_info()
            # Make a config object for each detector
            # In future check if a detector config changes during a run and 
            # create coords for each attr that changes in steps during the run.
            #axdat[det].coords['steps'] = range(nsteps)

            config_info = detector._xarray_info.get('attrs')
#            if attrs:
#                axdat[det].coords[det+'_config'] = ([det+'_steps'], range(nsteps))
#                axdat[det].coords[det+'_config'].attrs.update(attrs)

            for attr, attr_func in det_func.items():
                alias = attr_func.get('alias')
                attrs_info = attr_func.get('attr_info', {})
                if detector._tabclass == 'evtData':
                    if detector.evtData is not None:
                        item = detector.evtData._attr_info.get(attr)
                        if item:
                            attrs_info.update({a: item[a] for a in ['doc', 'unit']})
                    else:
                        print 'No data for {:} in {:}'.format(str(detector), attr)
                
                else:
                    if detector._calib_class is not None:
                        item = detector.calibData._attr_info.get(attr)
                        if item is not None:
                            attrs_info.update(item)
                    
                    elif detector.detector is not None:
                        item = detector.detector._attr_info.get(attr)
                        if item is not None:
                            attrs_info.update(item)

                if src_info:
                    attrs_info.update(src_info)

                # Make sure no None attrs
                for a, aitm in attr_info.items():
                    if aitm is None:
                        attr_info.update({a, ''})

                if 'event' in attr_func:
                    axdat[det][alias].attrs.update(attrs_info)
                
                if 'welford' in attr_func:
                    attr_func['welford']['attrs_info'] = attrs_info
                    attr_func['welford']['config_info'] = config_info

            coords = detector._xarray_info.get('coords')
            if coords:
                for coord, item in sorted(coords.items(), key=operator.itemgetter(0)):
                    try:
                        if isinstance(item, tuple):
                            dims = [det+'_'+dim for dim in item[0]]
                            vals = item[1]
                            axdat[det].coords[det+'_'+coord] = (dims,vals)
                        else:
                            axdat[det].coords[det+'_'+coord] = item
                    except:
                        print det, coord, item

    ds.reload()
    print 'xarray Dataset configured'

    time0 = time.time()
    igood = -1
    aievt = {}
    aievents = {}

    # keep track of events for each det
    for srcstr, srcitem in ds.configData._sources.items():
        det = srcitem.get('alias')
        aievt[det] = -1
        aievents[det] = []
   
    if ichunk is not None:
        ievent0 = (ichunk-1)*nevents
        print 'Making chunk {:}'.format(ichunk)
        print 'Starting with event {:} of {:}'.format(ievent0,ds.nevents)
        print 'Analyzing {:} events'.format(nevents)
        xbase.attrs['ichunk'] = ichunk
        if ichunk > 1:
            for i in range(ievent0):
                evt = ds.events.next()
        
            print 'Previous event before current chunk:', evt

    if ichunk is not None:
        evtformat = '{:10.1f} sec, Event {:} of {:} in chunk with {:} accepted'
    else:
        evtformat = '{:10.1f} sec, Event {:} of {:} with {:} accepted'
    
    #for ievent in range(ds.nevents+1):
    for ievt in range(nevents):
        ievent = ievent0+ievt
        if ievt > 0 and (ievt % 100) == 0:
            print evtformat.format(time.time()-time0, ievt, nevents, igood+1)
        
        if ievent < ds.nevents:
            try:
                evt = ds.events.next(publish=publish, init=publish)
            except:
                ievent = -1
                continue
        else:
            ievent = -1
            continue

        if len(set(eventCodes) & set(evt.Evr.eventCodes)) == 0:
            continue
       
        dtime = evt.EventId
        if dtime is None:
            continue
        
        igood += 1
        if igood+1 == nevents:
            break

        istep = ds._istep
        xbase['step'][igood] = istep
        btimes.append(dtime)
        for ec in evt.Evr.eventCodes:
            if ec in eventCodes:
                xbase['ec{:}'.format(ec)][igood] = True

        for attr, codes in code_flags.items():
            if evt.Evr.present(codes):
                xbase.coords[attr][igood] = True

        
        for det in evt._attrs:
            detector = evt._dets.get(det)
            atimes[det].append(dtime)
            aievt[det] += 1 
            ievt = aievt[det]
            aievents[det].append(ievent)
            
            for attr, attr_func in det_funcs.get(det, {}).items():
                vals = getattr(detector, attr)
                alias = attr_func.get('alias')
                if vals is not None and 'event' in attr_func:
                    # Fill event
                    try:
                        axdat[det][alias][ievt] = vals
                    except:
                        #print 'Event Error', alias, det, attr, ievent, vals
                        vals = None

                if vals is not None and 'welford' in attr_func:
                    # Add event to appropriate detector summary
                    for iec, fec in enumerate(attr_func['welford']['funcs'][istep]):
                        if eventCodes[iec] in evt.Evr.eventCodes:
                            try:
                                if not fec.shape or fec.shape == vals.shape: 
                                    fec(vals)
                                else:
                                    print 'Summary Add Error', alias, det, attr, \
                                        ievent, ec, fec.shape, vals.shape
                            except:
                                print 'Summary Error', alias, det, attr, \
                                        ievent, ec, fec, vals

    xbase = xbase.isel(time=range(len(btimes)))
    xbase['time'] =  [e.datetime64 for e in btimes]
    for attr, dtyp in ttypes.items():
        xbase.coords[attr] = (['time'], np.array([getattr(e, attr) for e in btimes],dtype=dtyp))
        
    # fill each control PV with current step value
    scan_variables = []
    if ds.configData.ScanData and ds.configData.ScanData.nsteps > 1:
        for attr, vals in ds.configData.ScanData.control_values.items():
            alias = ds.configData.ScanData.pvAliases[attr]
            scan_variables.append(alias)
            xbase.coords[alias] = (['time'], xbase.coords[alias+'_steps'][xbase.step]) 

    xbase.attrs['scan_variables'] = scan_variables
    xbase.attrs['correlation_variables'] = []

    #return xbase, axdat

    # cut down size of xdat
    det_list = [det for det in axdat]
    for det in sort(det_list):
        nevents = len(atimes[det])
        if nevents > 0 and det in axdat:
            #try:
                print 'merging', det
                xdat = axdat.pop(det)
                if 'time' in xdat:
                    xdat = xdat.isel(time=range(nevents))
                    xdat['time'] = [e.datetime64 for e in atimes[det]]
                    xdat = xdat.reindex_like(xbase)
        
                _make_summary(xdat, {det: det_funcs[det]})
                if make_images:
                    _make_images(xdat, {det: det_funcs[det]})
                
                xbase = xbase.merge(xdat)
            #except:
            #    print 'Could not merge', det

    attrs = [attr for attr,item in xbase.data_vars.items()] 
    for attr in attrs:
        for a in ['unit', 'doc']:
            if a in xbase[attr].attrs and xbase[attr].attrs[a] is None:
                xbase[attr].attrs[a] = ''

    return xbase
 
def _make_summary(xbase, det_funcs):
    for det, det_func in det_funcs.items():
        for attr, attr_func in det_func.items():
            alias = attr_func.get('alias')
            if 'welford' in attr_func:
                attrs_info = attr_func['welford']['attrs_info']
                nsteps = attr_func['welford']['shape'][0]
                neventCodes = attr_func['welford']['shape'][1]
                aevents = np.zeros(shape=(nsteps,neventCodes))
                for istep in range(nsteps):
                    for iec, fec in enumerate(attr_func['welford']['funcs'][istep]):
                        aevents[istep, iec] = fec.n

                print 'Events', det, attr_func['welford']['dims'][0:2], aevents
                xbase[det+'_events'] = (attr_func['welford']['dims'][0:2], aevents)
                xbase[det+'_events'].attrs.update(attr_func['welford']['config_info'])

                sum_dims = attr_func['welford']['dims']
                for fattr in attr_func['welford']['attrs']:
                    asums = np.zeros(shape=attr_func['welford']['shape'])
                    for istep in range(nsteps):
                        for iec, fec in enumerate(attr_func['welford']['funcs'][istep]):
                            vals = getattr(fec, fattr)()
                            if np.any(vals):
                                asums[istep, iec] = vals         
                    
                    print 'summary', det, fattr, sum_dims, np.array(asums).shape
                    xbase[alias+'_'+fattr] = (sum_dims, asums)
                    xbase[alias+'_'+fattr].attrs.update(attrs_info)
                
def _make_images(xbase, det_funcs):
    for det, det_func in det_funcs.items():
        for attr, attr_func in det_func.items():
            alias = attr_func.get('alias')
            if 'welford' in attr_func:
                attrs_info = attr_func['welford']['attrs_info']
                # Add event to appropriate detector summary
                if hasattr(xbase, det+'_indexes_x') and hasattr(xbase, det+'_indexes_y'):
                    indexes_x = np.array(getattr(xbase, det+'_indexes_x'))
                    indexes_y = np.array(getattr(xbase, det+'_indexes_y'))
                    img_dims = ['steps', 'codes', det+'_ximage', det+'_yimage']
                else:
                    continue

                print det, attr, alias
                try:
                    for fattr in attr_func['welford']['attrs']:
                        aimgs = []
                        nsteps = len(attr_func['welford']['funcs'])
                        for istep in range(nsteps):
                            imgs = []
                            for iec, fec in enumerate(attr_func['welford']['funcs'][istep]):
                                vals = getattr(fec, fattr)()
                                if np.any(vals):
                                    img = map_indexes(indexes_x, indexes_y, vals)
                                    #img = img_from_pixel_arrays(indexes_x, indexes_y, vals)
                                else:
                                    img = np.zeros((indexes_x.max()+1,indexes_y.max()+1))
                                
                                imgs.append(np.array(img))
                        
                            aimgs.append(imgs)
                        
                        print 'image', det, fattr, img_dims, np.array(aimgs).shape
                        xbase[det+'_image_'+fattr] = (img_dims, np.array(aimgs))
                        xbase[det+'_image_'+fattr].attrs.update(attrs_info)
                except:
                    print 'No images possible for', attr, alias

# simplified from img_from_pixel_arrays
def map_indexes(xx, yy, ww):                                                                      
    a = np.zeros([xx.max()+1,yy.max()+1])
    a[xx,yy] = ww
    return a

#class Dataset(xarray.Dataset):
#
#    def init(self, ds=None, **kwargs):
#
#        xarray.Dataset.__init__()
#
#    def make_foo(self):
#        self.attrs['foo'] = 'bar'
#


