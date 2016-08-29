import time
import os
import logging

import pandas as pd
from pylab import *
import xarray as xr
import seaborn as sns

import output_html
import markup

pd.set_option('precision',2)
pd.set_option('max_colwidth', 66)
pd.set_option('display.width', 110)
plt.rcParams['axes.labelsize'] = 16 

def hover(hover_color="#ffff99"):
    return dict(selector="tr:hover",
                props=[("background-color", "%s" % hover_color)])

styles = [
    hover(),
    dict(selector="th", props=[("font-size", "100%"),
                               ("text-align", "center")]),
    dict(selector="caption", props=[("caption-side", "bottom")])
]

def get_run_size(run, start_path = '.'):
    tt = 0
    nn = 0
    runstr = 'r{:04.0f}'.format(run)
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            if 'xtc' in f.lower() and runstr in f.lower():
                fp = os.path.join(dirpath,f)
                tt += os.path.getsize(fp)
                nn += 1
    return nn, tt

def load_DataSource(**kwargs):
    import PyDataSource
    return PyDataSource.DataSource(**kwargs)

def load_xdat(ds, **kwargs):
    from psxarray import get_xdat
    get_xdat(ds, **kwargs)

# Make tidy pandas dataframe
#a = x.dims.keys()
#a.remove('time')
#df = x.drop(a).to_dataframe()

class Build_html(object):
    """Class to build an html RunSummary report.
    """

    def __init__(self, xdat=None, ds=None, **kwargs):
        
        self.logger = logging.getLogger(__name__+'.build_html')
        self.logger.info(__name__)
       
        self.results = {}

        if xdat:
            if xdat.__class__.__name__ == 'DataSource':
                # 1st arg is actually PyDataSource.DataSource
                ds = xdat
                xdat = None

        elif not ds:
            if 'run' in kwargs and 'exp' in kwargs:
                self.logger.info('Loading PyDataSource data')
                ds = load_DataSource(**kwargs)
            else:
                print 'Require valid xdat=xarray object and/or ds=PyDataSource object'
                print 'or alternatively exp and run kwargs'
                return

        if xdat:
            self._xdat = xdat
        else:
            self.logger.info('Building xarray data')
            self._xdat = load_xdat(ds, **kwargs)
           
        self.exp = str(self._xdat.experiment)
        self.instrument = str(self._xdat.instrument)
        self.run = int(self._xdat.run.values[0])
        if ds:
            self._ds = ds
        else:
            self.logger.info('Loading PyDataSource data')
            self._ds = load_DataSource(exp=self.exp, run=self.run)

        self.title = title='{:} Run {:}'.format(self.exp ,self.run)
       
        self.xtc_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp,'xtc')
        self.res_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp,'res')
        self.scratch_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp,'scratch')
        self.nc_file = '{:}/nc/run{:04}.nc'.format(self.scratch_dir, self.run)
       
        self.aliases = list(set([item.attrs.get('alias') for attr,item in self._xdat.data_vars.items() if item.attrs.get('alias')]))

        if 'steps' in self._xdat:
            self.nsteps = len(self._xdat.steps)
        else:
            self.nsteps = None

        if 'codes' in self._xdat:
            self.ncodes = len(self._xdat.codes)
        else:
            self.ncodes = None

        self._init_output()
    
    def _init_output(self, path=None, **kwargs):
        """Set output path and build appropriate directories for html
        """
        self._nfile, self._nbytes = get_run_size(self.run, start_path=self.xtc_dir)
        self.logger.info( "counting files in {:}".format(self.xtc_dir) )
        
        if not path:
            path = os.path.join(self.res_dir, 'RunSummary')
        elif path == 'home':
            path = os.path.join(os.path.expanduser('~'), 'RunSummary', self.exp)
        
        print 'Setting output path to', path

        self.path = path
        self.output_dir = os.path.join(self.path, 'run{:04}'.format(self.run))
       
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
 
  
    def add_all(self, quiet=True, **kwargs):
        """Add All detectors.
        """
        # Add RunStats
        self.add_detstats()

        for alias in self.aliases:
            if not quiet:
                print 'adding detector', alias
            if 'cuts' in self._xdat.attrs:
                for cut in self._xdat.attrs['cuts']:
                    self.add_detector(alias, cut=cut, **kwargs)
                
            else:
                self.add_detector(alias, **kwargs)

        plt.cla()

    def add_detstats(self, alias='RunStats'):
        """Add detector statistics. 
        """
        attrs = [a for a in self._xdat.variables.keys() if a.endswith('events')]
        if attrs and self.nsteps == 1:
            if alias not in self.results: 
                self.results[alias] = {'figure': {}, 'table': {}, 'text': {}}
        
            df_tbl = self._xdat[attrs].sel(steps=0).to_array().to_pandas()
            self.eventStats = df_tbl
            doc = []
            doc.append("A summary of the mean, stdev, min and max values were created "
                      +"for each eventCode for the detectors in this table:")
            howto = []
            howto.append("attrs={:})".format(attrs))
            howto.append("df_tbl = x[attrs].sel(steps=0).to_array().to_pandas()")
            self.results[alias]['table'].update({alias:{'DataFrame': df_tbl, 
                                                        'name': 'df_tbl',
                                                        'howto': howto, 
                                                        'doc': doc}})

    def add_detector(self, alias='FEEGasDetEnergy', percentiles=[0.05,0.50,0.95], 
                           attrs=None,
                           figsize=None, 
                           layout=None,
                           groupby=None,
                           scat_name=None,
                           make_scatter=False,
                           make_table=True,
                           make_timeplot=True,
                           make_histplot=True,
                           labelsize=20,
                           robust_attrs=None,
                           plot_errors=False,
                           cut=None, 
                           show=False):
        """Add a single detector based on alias
        """
        x = self._xdat
        #if 'Damage_cut' in x:
        #    x = x.where(x.Damage_cut == 1)

        plt.rcParams['axes.labelsize'] = labelsize 
        desc_attrs = ['unit','doc']
        default_scatter_attrs = {
                'EBeam': ['ebeamCharge', 'ebeamDumpCharge','ebeamL3Energy', 'ebeamPhotonEnergy', 'ebeamUndPosX', 'ebeamUndPosY']
                #'EBeam': ['ebeamCharge', 'ebeamL3Energy', 'ebeamPhotonEnergy', 'ebeamXTCAVAmpl', 'ebeamXTCAVPhase']
                }

        if groupby is None:
            groupby = x.attrs.get('scan_variables')

        if groupby and not np.shape(groupby):
            groupby = [groupby]

        attr_names = None
        if not attrs:
            attrs = [attr for attr,item in x.data_vars.items() \
                     if item.attrs.get('alias') == alias and item.dims == ('time',)]
            attrs0 = [attr for attr in attrs]
            for attr in x.attrs.get('correlation_variables', []):
                attrs.append(attr)
        else:
            if isinstance(attrs, dict):
                attr_names = attrs
                attrs = attr_names.values()
            
            attrs0 = [attr for attr in attrs]
        
        if groupby:
            for group in groupby:
                attrs.append(group)

        if not attr_names:
            attr_names = {attr: attr.replace(alias+'_','') for attr in attrs if attr in x}

        catagory = alias
#        if cut is None:
#            if 'Damage_cut' not in x:
#                self.make_damage_cut()
#            cut = 'Damage_cut'
        
        if cut:
            catagory = alias+' '+cut

        if attrs:
            if catagory not in self.results: 
                self.results[catagory] = {'figure': {}, 'table': {}, 'text': {}}
            
            nattrs = len(attrs)
            if nattrs == 4:
                nrows = 2
                ncolumns = 2
            elif nattrs == 2:
                nrows = 2
                ncolumns = 1
            else:
                ncolumns = int(min([nattrs,3]))
                nrows = int(np.ceil(nattrs/float(ncolumns)))

            if not layout:
                layout = (nrows,ncolumns)
            if not figsize:
                figsize = (ncolumns*4,nrows*3.5)
            
            print attrs
            xselect = x[attrs]
            tselect = 'time'
            if cut:
                tselect = 'points'
                itimes = xselect.groupby(cut).groups
                if len(itimes) > 1:
                    xselect = xselect.isel_points(time=itimes[1])

            df = xselect.to_array().to_pandas().T
            
            howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
            howto.append("attrs = {:}".format(attrs))
            if cut:
                howto.append("itimes = x[attrs].groupby({:}).groups[1]".format(cut))
                howto.append("xselect = x[attrs].isel_points(times=itimes)")
            else:
                howto.append("xselect = x[attrs]")

            howto.append("df = xselect.to_array().to_pandas().T")

            df_attrs = pd.DataFrame({attr_names[attr]: {a: x[attr].attrs.get(a) for a in desc_attrs} \
                                     for attr in attrs}).T[desc_attrs]
            
            # Need to add in howto print formatted df_attrs, but tricky
            self.results[catagory]['table'].update({'attrs': {'DataFrame': df_attrs, 'name': 'df_attrs',
                        'format': {'unit': '{:<10s}', 'doc': '{:<69s}'}}})
            
            df.rename(inplace=True, columns=attr_names)
            howto.append("attr_names={:}".format(attr_names))
            howto.append("df.rename(inplace=True, columns=attr_names)")
            self.results[catagory]['text'].update({'setup':{'howto': howto}})
            
            #df_attrs.to_string(formatters={'unit': '{:<10s}'.format, 'doc': '{:<60s}'.format},justify='left')

            df_tbl = df.describe(percentiles=percentiles).T.round({'count':0})
            howto = ["df_tbl = df.describe(percentiles={:}).T.round({:})".format(percentiles,{'count':0})]
            howto.append("print df_tbl")

            # make table using pandas
            if make_table:
                self.results[catagory]['table'].update({'stats':{'DataFrame': df_tbl, 'name': 'df_tbl', 'howto': howto}})
           
            # make time plots
            if make_timeplot:
                plt_type = 'time'
                df.plot(subplots=True, sharex=True, layout=layout, figsize=figsize)
                howto = ["df.plot(subplots=True, sharex=True, layout={:}, figsize={:})".format(layout, figsize)]
                self.add_plot(catagory, plt_type, howto=howto)

            # Make groupby plots
            if groupby:
                if not isinstance(groupby, list):
                    groupby = [groupby]
                for grp in groupby:
                    
#                    try:
                    if True:
                        # not sure why this is a list sometimes
                        if isinstance(grp, list):
                            grp = grp[0]

                        group = str(grp)
                        print group, tselect
                        if plot_errors:
                            df_group = x[attrs].to_array().to_pandas().T.groupby(group)
                            xaxis = df_group[group].mean().values
                            xlab = group
                            unit = x[group].attrs.get('unit')
                            if unit:
                                xlab = '{:} [{:}]'.format(xlab, unit)
                            
                            for attr in attrs0:
                                plt_type = '{:} vs {:}'.format(attr, group)
                                plt.figure()
                                yaxis = df_group[attr].mean().values
                                yerr = df_group[attr].std().values
                                plt.errorbar(xaxis,yaxis,yerr=yerr)
                                plt.xlabel(xlab)
                                ylab = attr
                                unit = x[attr].attrs.get('unit')
                                if unit:
                                    ylab = '{:} [{:}]'.format(ylab, unit)
                                
                                plt.ylabel(ylab)
                                howto = ["xaxis = x['{:}'].values".format(group)]
                                howto.append("df_group = xselect.to_array().to_pandas().T.groupby('{:}')".format(group))
                                howto.append("xaxis = df_group['{:}'].mean()".format(group)) 
                                howto.append("yaxis = df_group['{:}'].mean()".format(attr)) 
                                howto.append("yerr = df_group['{:}'].std()".format(attr)) 
                                howto.append("figure()")
                                howto.append("plt.errorbar(xaxis,yaxis,yerr=yerr)")
                                howto.append("plt.xlabel('{:}')".format(xlab))
                                howto.append("plt.ylabel('{:}')".format(ylab))
                                self.add_plot(catagory, plt_type, howto=howto)

                        else: 
                            plt_type = 'vs {:}'.format(group)
                            x_group = xselect.groupby(group).mean(dim=tselect)
                            df_group = x_group.to_array().to_pandas().T
                            df_group.plot(subplots=True, sharex=True, layout=layout, figsize=figsize)
                            howto = ["x_group = xselect.groupby('{:}').mean(dim={:})".format(group,tselect)]
                            howto.append("df_group = x_group.to_array().to_pandas().T")
                            howto.append("df_group.plot(subplots=True, sharex=True, layout={:}, figsize={:})".format(layout, figsize))
                            self.add_plot(catagory, plt_type, howto=howto)
                        
                        #grp_values = set(xselect[grp].values)
                        #if len(grp_values) < 9:
                        #    for attr in attrs0:
                        #        plt.cla()

                        #        xselect[attr] plot(robust=True)
                        #        howto = []
                        #        howto.append("xdat.squeeze('steps').sel(codes={:}).plot(robust=True)".format(code))
                        #        howto.append("plt.show()")
                        #        name = attr.replace(alias+'_','')
                        #        plt_type = 'summary_{:}_ec{:}'.format(name, code) 
                        #        self.add_plot(alias, plt_type, howto=howto)

                    else:
                    #except:
                        print 'Failed adding:', group
                        plt.cla()
                        plt.close()

            # make hist plots:
            if make_histplot:
                #perlow = '{:}%'.format(int(min(percentiles)*100))
                #perhigh = '{:}%'.format(int(max(percentiles)*100))
                #dfcut = df[(df > df_tbl[perlow]-2*df_tbl['std']).all(axis=1) & (df < df_tbl[perhigh]+2*df_tbl['std']).all(axis=1)]
                dfcut = df[(df > df_tbl['5%']-2*df_tbl['std']).all(axis=1) & (df < df_tbl['95%']+2*df_tbl['std']).all(axis=1)]
                try:
                    howto = ["dfcut = df[(df > df_tbl['5%']-2*df_tbl['std']).all(axis=1) & (df < df_tbl['95%']+2*df_tbl['std']).all(axis=1)]"] 
                    plt_type = 'hist'
                    dfcut.hist(alpha=0.2, layout=layout, figsize=figsize)
                    howto.append("dfcut.hist(alpha=0.2, layout={:}, figsize={:})".format(layout, figsize))
                    self.add_plot(catagory, plt_type, howto=howto)
                except:
                    plt.cla()
                    plt.close()
                    print 'make_histplot failed'
                    print howto
                    print dfcut.keys()

            # make scatter matrix using pandas -- cut outliers
            if make_scatter:
                if not robust_attrs:
                    robust_attrs = df.keys()
                    
                print df_tbl.keys()
                print df_tbl.T.keys()
                print df.keys()
                print robust_attrs
                #return df, robust_attrs, df_tbl
                robust_attrs = [a for a in [attr.replace(alias+'_','') for attr in robust_attrs] if a in df.keys()]
                dfr = df[robust_attrs]
                df_tblr = df_tbl.T[robust_attrs].T
                dfcut = df[(dfr > df_tblr['5%']-2*df_tblr['std']).all(axis=1) & (dfr < df_tblr['95%']+2*df_tblr['std']).all(axis=1)]

                if isinstance(make_scatter, list):
                    #scat_attrs = [attr.replace(alias+'_','') for attr in make_scatter if attr in dfcut.keys()]
                    scat_attrs = [a for a in [attr.replace(alias+'_','') for attr in make_scatter] if a in dfcut.keys()]

                else:
                    scat_attrs = default_scatter_attrs.get(alias, attr_names.values())
                
                for attr in x.attrs.get('correlation_variables'):
                    if attr in dfcut.keys() and attr not in scat_attrs:
                        scat_attrs.append(attr)

                try:
                    print scat_attrs
                    print dfcut.keys()
                    dfscat = dfcut[scat_attrs]
                except:
                    plt.cla()
                    plt.close()
                    print 'make_scatter failed'
                    return dfcut

                try:
                    howto = ["scat_attrs = {:}".format(scat_attrs)]
                    howto.append("dfscat = dfcut[scat_attrs]")
                    if groupby:
                        if scat_name:
                            plt_type = scat_name
                        else:
                            plt_type = 'correlation with {:}'.format(group)
                        
                        pltattrs = scat_attrs
                        for group in groupby:
                            if group in pltattrs:
                                pltattrs.remove(group)
                            
                        sns.set()
                        plt.rcParams['axes.labelsize'] = 10 
                        g = sns.pairplot(dfcut, hue=group,
                                x_vars=pltattrs,y_vars=pltattrs,
                                #palette="Set1", 
                                size=2.5) 
                                #palette="Set1", diag_kind="kde", size=2.5) 
                                #x_vars=pltattrs,y_vars=pltattrs)
                        #g = sns.PairGrid(dfscat, hue=group, 
                        #        palette="Set1", 
                        #        x_vars=pltattrs,y_vars=pltattrs)
                        #g = g.map_offdiag(plt.scatter)
                        #g = g.add_legend()
                        self.add_plot(catagory, plt_type, howto=howto)
#df = x[attrs].to_dataframe()
#sns.pairplot(df, hue="LaserOn", vars=attrs0)
                    else:

                        plt_type = 'scatter_matrix'
                        Axes = pd.tools.plotting.scatter_matrix(dfscat, alpha=0.2, figsize=(20, 20), diagonal='kde')
                        howto.append("Axes = pd.tools.plotting.scatter_matrix(dfscat, alpha=0.2, figsize=(20, 20), diagonal='kde')")
                        self.add_plot(catagory, plt_type, howto=howto)
    
    #                    #y ticklabels
    #                    [plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
    #                    howto.append("[plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]")
    #                    #x ticklabels
    #                    [plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
    #                    howto.append("[plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]")
    #                    #y labels
    #                    [plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]
    #                    howto.append("[plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]")
    #                    #x labels
    #                    [plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]
    #                    howto.append("[plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]")

                except:
                    plt.cla()
                    plt.close()
                    print dfscat

        plt.close('all')
   
    def add_scatter(self, df, catagory='scatter', doc=None, group=None, attrs=None, howto=[]):
        if not attrs:
            attrs = df.keys()
       
        sns.set()
        plt.rcParams['axes.labelsize'] = 10 
        pltattrs = [attr for attr in attrs if attr not in [group]]

        howto.append("sns.set()")
        howto.append("pltattrs = {:}".format(pltattrs))
        
        if group:
            plt_type = 'correlation with {:}'.format(group)
            g = sns.pairplot(df, hue=group,
                    x_vars=pltattrs,y_vars=pltattrs,
                    size=2.5) 
            howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")
        else:
            plt_type = 'correlation'
            g = sns.PairGrid(df, x_vars=pltattrs,y_vars=pltattrs, size=2.5) 
            g = g.map_upper(plt.scatter)
            g = g.map_lower(sns.kdeplot, cmap="Blues_d")
            g = g.map_diag(sns.kdeplot, lw=3, legend=False)

            howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")

        self.add_plot(catagory, plt_type, howto=howto)
        plt.close('all')

    def add_scatter_groups(self, attr_groups=None, group=None, howto=None, attrs=None):
        """
        """
        x = self._xdat
        if not attrs:
            attrs=['EBeam_ebeamPhotonEnergy', 'FEEGasDetEnergy_f_21_ENRC']
            for attr in x.correlation_variables:
                attrs.append(attr)

        if not group:
            group = 'XrayOn'
        
        if group not in x:
            print group, 'is not a valid group -- ignoring group keyword'
            group = None

        if not howto:
            howto = []
        
        if not attr_groups:
            attr_groups = {'Undulator X': ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX'],
                           'Undulator Y': ['EBeam_ebeamUndAngY', 'EBeam_ebeamUndPosY'],
                           'LTU X': ['EBeam_ebeamLTUAngX', 'EBeam_ebeamLTUPosX'],
                           'LTU Y': ['EBeam_ebeamLTUAngY', 'EBeam_ebeamLTUPosY'],
                           'LTU 250 and 450': ['EBeam_ebeamLTU250', 'EBeam_ebeamLTU450'],
                           'Energy BC': ['EBeam_ebeamEnergyBC1', 'EBeam_ebeamEnergyBC2'],
                           'Charge': ['EBeam_ebeamCharge', 'EBeam_ebeamDumpCharge'],
                           'Peak Current': ['EBeam_ebeamPkCurrBC1', 'EBeam_ebeamPkCurrBC2'],
                           'XTCAV': ['EBeam_ebeamXTCAVAmpl', 'EBeam_ebeamXTCAVPhase'],
                           'Phase Cavity Charge': ['PhaseCavity_charge1', 'PhaseCavity_charge2'],
                           'Phase Cavity Fit': ['PhaseCavity_fitTime1', 'PhaseCavity_fitTime2'],
                           'Gasdet': ['FEEGasDetEnergy_f_63_ENRC', 'FEEGasDetEnergy_f_11_ENRC'],
                           }

        all_attrs = [attr for attr,item in x.data_vars.items() if item.dims == ('time',)]
        
        if 'Damage_cut' not in x:
            self.make_damage_cut()

        df = self._xdat[all_attrs].where(self._xdat.Damage_cut == 1).to_array().to_pandas().T
        
        if group not in df:
            df[group] = x[group].values

        df = df.dropna()

        for catagory, grp_attrs in attr_groups.items():
            gattrs = [attr for attr in attrs]
            for attr in grp_attrs:
                gattrs.append(attr)
            ghowto = howto
            print catagory, gattrs, group, ghowto
            self.add_scatter(df, catagory=catagory, attrs=gattrs, howto=ghowto, group=group)

    def make_damage_cut(self, charge_min=1., 
            fitTime1_min=0.3, fitTime1_max=1.1, 
            fitTime2_min=-0.4, fitTime2_max=0.2,
            gasdet_min=0.5):
        """
        """
        x = self._xdat
        phasecut = (x.PhaseCavity_charge1.values > charge_min) & \
                (x.PhaseCavity_charge2.values > charge_min) & \
                (x.PhaseCavity_fitTime1.values <= fitTime1_max) & \
                (x.PhaseCavity_fitTime1.values >= fitTime1_min) & \
                (x.PhaseCavity_fitTime2.values <= fitTime2_max) & \
                (x.PhaseCavity_fitTime2.values >= fitTime2_min) 
        gasdetcut =  x.FEEGasDetEnergy_f_11_ENRC.values > 0.5
        damagecut = phasecut & gasdetcut & (x.EBeam_damageMask.values == 0)  
        self._xdat.coords['Damage_cut'] = (['time'], damagecut)

#    def add_scatter(self, df, catagory='scatter', group=None, attrs=None, howto=[]):
#        if not attrs:
#            attrs = df.keys()
#       
#        pltattrs = [attr for attr in attrs if attr not in [group]]
#        if group:
#            plt_type = 'correlation with {:}'.format(group)
#        else:
#            plt_type = 'correlation'
#
#        sns.set()
#        plt.rcParams['axes.labelsize'] = 10 
#        g = sns.pairplot(df, hue=group,
#                x_vars=pltattrs,y_vars=pltattrs,
#                size=2.5) 
#
#        howto.append("sns.set()")
#        howto.append("pltattrs = {:}".format(pltattrs))
#        howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")
#
#        self.add_plot(catagory, plt_type, howto=howto)
#
#                #palette="Set1", diag_kind="kde", size=2.5) 
#                #x_vars=pltattrs,y_vars=pltattrs)
#        #g = sns.PairGrid(dfscat, hue=group, 
#        #        palette="Set1", 
#        #        x_vars=pltattrs,y_vars=pltattrs)
#        #g = g.map_offdiag(plt.scatter)
#        #g = g.add_legend()
 

    def add_plot(self, catagory, plt_type, howto=[], show=False):
        if catagory not in self.results:
            self.results[catagory] = {'figure': {}, 'table': {}, 'text': {}}
        
        plt_file = '{:}_{:}.png'.format(catagory, plt_type).replace(' ','_') 
        self.results[catagory]['figure'].update({plt_type: {'path': self.output_dir, 
                                                         'png': plt_file,
                                                         'howto': howto}})
        plt.savefig(os.path.join(self.output_dir, plt_file))
        if show:
            plt.show()
        else:
            plt.close()

    def add_summary(self, 
                        variables=None, 
                        max_columns=10,
                        show=False,
                        groupby=None,
                        codes=None,
                        wrap=False,
                        labelsize=None, 
                        figsize=None,
                        layout=None):

        if not variables:
            if 'scan_variables' in self._xdat.attrs:
                variables = self._xdat.scan_variables
            else:
                variables = []

        plt.close()
        plt.cla()
        variables = [v for v in variables if v in self._xdat.variables.keys()]

        if groupby is None:
            groupby = self._xdat.attrs.get('scan_variables')

        if groupby and not variables:
            variables = [attr for attr, item in self._xdat.variables.items() \
                            if len(item.shape) == 2 and item.shape[1] > 8]

        for attr in variables:
            if 'alias' in self._xdat[attr].attrs:
                alias = self._xdat[attr].alias
            else:
                alias = attr.split('_')[0]

            if alias not in self.results:
                self.results[alias] = {'figure': {}, 'table': {}, 'text': {}}

            # Make groupby plots
            if groupby:
                if not isinstance(groupby, list):
                    groupby = [groupby]
                for group in groupby:
                    xdat = self._xdat[attr].T
                    x_group = xdat.groupby(group).mean(dim='time')
                    if len(x_group.shape) == 2:
                        x_group.plot()
                        plt_type = '{:} vs {:}'.format(attr.lstrip(alias+'_'), group) 
                        howto = ["x_group = x['{:}'].T.groupby('{:}').mean(dim='time')".format(attr, group)]
                        howto.append('x_group.plot()')
                        self.add_plot(alias, plt_type, howto=howto)

            if len(self._xdat[attr].shape) == 2:
#                self._xdat[attr].mean(axis=0).plot()
#                plt_type = '{:} summary'.format(attr.lstrip(alias+'_')) 
#                howto = ["x['{:}'].mean(axis=0).plot()".format(attr)]
#                self.add_plot(alias, plt_type, howto=howto)
                
                self._xdat[attr].plot()
                plt_type = '{:} with {:}'.format(attr.lstrip(alias+'_'), 'time') 
                howto = ["x['{:}'].plot()".format(attr)]
                self.add_plot(alias, plt_type, howto=howto)

            if len(self._xdat[attr].shape) in [2,3]:
                for iaxis in range(3):
                    pltaxis = self._xdat[attr].dims[iaxis]
                    self._xdat[attr].sum(axis=iaxis).plot()
                    plt_type = '{:} sum over {:}'.format(attr.lstrip(alias+'_'), pltaxis) 
                    howto = ["x['{:}'].mean(axis=0).plot()".format(attr)]
                    self.add_plot(alias, plt_type, howto=howto)
                

            if len(self._xdat[attr].shape) == 4 and 'codes' in self._xdat[attr]:
                xdat = self._xdat[attr]
                attr_codes = None
                if codes:
                    attr_codes = codes
                elif hasattr(self, 'eventStats'):
                    try:
                        tbl = self.eventStats.T.get(alias+'_events')
                        if tbl is not None:
                            attr_codes = [c for c,n in tbl.to_dict().items() if n > 0]
                    except:
                        print 'WARNING: Need to fix auto getting attr_codes'
                
                if not attr_codes:
                    if 'codes' in self._xdat:
                        attr_codes = self._xdat.codes.values
                    else:
                        attr_codes = [140]

                if self.nsteps < max_columns**2:
                    ncolumns = int(sqrt(self.nsteps))
                else:
                    ncolumns = int(min([self.nsteps,max_columns]))
                    
                nrows = int(np.ceil(self.nsteps/float(ncolumns)))
                
                if not labelsize:
                    if ncolumns > 9:
                        labelsize = 12
                    elif ncolumns > 5:
                        labelsize = 18
                    else:
                        labelsize = 24

                if not layout:
                    layout = (nrows,ncolumns)
                if not figsize:
                    figsize = (ncolumns*4,nrows*4)
                
                plt.rcParams['axes.labelsize'] = 24 
                howto = ["plt.rcParams['axes.labelsize'] = {:}".format(24)]
                howto.append("attr = '{:}'".format(attr))
                howto.append("xdat = x[attr]")
                self.results[alias]['text'].update({'setup':{'howto': howto}})

                # Need to fix errors calculation with std and plot with error bars
                #fig, ax = plt.subplots()
                #errors = sqrt(self._xdat[attr].sum(axis=(2,3)).to_pandas()**2)
                #df.plot(yerr=errors, ax=ax)

                if self.ncodes > 4:
                    col_wrap = int(min([self.ncodes,3]))
                    aspect = 0.9
                elif self.ncodes == 4:
                    col_wrap = 2
                    aspect = 0.9
                else:
                    col_wrap = 1
                    aspect = 0.9
                
                if self.nsteps == 1:
                    if wrap:
                        howto = []
                        df = xdat.squeeze('steps')
                        howto.append("df = xdat.squeeze('steps')")
                        if self.ncodes == 1:
                            df.squeeze('codes').plot(robust=True)
                            howto.append("df.squeeze('codes').plot()")
                            howto.append("plt.show()")
                        else:
                            df.plot(col='codes', col_wrap=col_wrap, size=20, aspect=aspect, robust=True)
                            howto.append("df.plot(col='codes', col_wrap={:})".format(col_wrap))
                            howto.append("plt.show()")

                        name = attr.replace(alias+'_','')
                        plt_type = 'summary_{:}'.format(name) 
                        self.add_plot(alias, plt_type, howto=howto)
                        
                    else:
                        for code in attr_codes:
                            plt.cla()
                            xdat.squeeze('steps').sel(codes=code).plot(robust=True)
                            howto = []
                            howto.append("xdat.squeeze('steps').sel(codes={:}).plot(robust=True)".format(code))
                            howto.append("plt.show()")
                            name = attr.replace(alias+'_','')
                            plt_type = 'summary_{:}_ec{:}'.format(name, code) 
                            self.add_plot(alias, plt_type, howto=howto)

                else:
                    # plot mean of image vs steps for each event code
                    df = xdat.sum(axis=(2,3)).to_pandas()
                    df.plot(subplots=True)
                   
                    howto = []
                    howto.append("df = xdat.sum(axis=(2,3)).to_pandas()")
                    howto.append("df.plot(subplots=True)")
                    howto.append("plt.show()")
                    name = attr.replace(alias+'_','')
                    plt_type = 'scan_{:}_summary'.format(name) 
                    self.add_plot(alias, plt_type, howto=howto)
     
                    for iaxis in [2,3]:
                        df = xdat.mean(axis=iaxis)
                        
                        howto = ["plt.rcParams['axes.labelsize'] = 24"]
                        howto.append("df = xdat.mean(axis={:}).to_pandas()".format(iaxis))
                        if self.ncodes == 1:
                            df.squeeze('codes').plot(robust=True)
                            howto.append("df.squeeze('codes').plot(robust=True)")
                            howto.append("plt.show()")
                        else:
                            df.plot(col='codes', col_wrap=col_wrap, size=20, aspect=aspect, robust=True)
                            sformat = "df.plot(col='codes', col_wrap={:}, size=20, aspect={:}, robust=True)"
                            howto.append(sformat.format(col_wrap, aspect))
                            howto.append("plt.show()")

                        name = df.dims[2].replace(alias+'_','')
                        plt_type = 'scan_vs_{:}'.format(name) 
                        self.add_plot(alias, plt_type, howto=howto)

                    plt.rcParams['axes.labelsize'] = labelsize
                    for code in self._xdat[attr].codes.values:
                        df = xdat.sel(codes=code)
                        df.plot(col='steps', col_wrap=ncolumns, robust=True)
                        
                        howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
                        howto.append("df = xdat.sel(codes={:})".format(code))
                        howto.append("df.plot(col='steps', col_wrap={:})".format(ncolumns))
                        name = attr.replace(alias+'_','')
                        plt_type = 'scan_{:}_2D_ec{:}'.format(name, str(code)) 
                        self.add_plot(alias, plt_type, howto=howto)

    def add_event(self, max_columns=10,
                        nmax_events=100,
                        show=False,
                        variables=None,
                        labelsize=None, 
                        figsize=None,
                        layout=None):
        """Add every event of variables.
        """

        if not variables:
            if 'event_variables' in self._xdat.attrs:
                variables = self._xdat.event_variables
            else:
                variables = []

        variables = [v for v in variables if v in self._xdat.variables.keys()]

        for attr in variables:
            alias = self._xdat[attr].alias
            if alias not in self.results:
                self.results[alias] = {'figure': {}, 'table': {}, 'text': {}}

            if len(self._xdat[attr].shape) == 3 and 'time' in self._xdat[attr]:
                xdat = self._xdat[attr]
                itimes = xdat.groupby('ec{:}'.format(xdat.eventCode)).groups[1]
                df = xdat.isel_points(time=itimes)
                #cut = [i for i,a in enumerate(self._xdat[attr][:,0,0].values) if isfinite(a)]
               
                nevents = len(df.time)
                if nevents < max_columns**2:
                    ncolumns = int(sqrt(nevents))
                else:
                    ncolumns = int(min([nevents,max_columns]))
                    
                nrows = int(np.ceil(nevents/float(ncolumns)))
                
                if not labelsize:
                    if ncolumns > 9:
                        labelsize = 12
                    elif ncolumns > 5:
                        labelsize = 18
                    else:
                        labelsize = 24

                if not layout:
                    layout = (nrows,ncolumns)
                if not figsize:
                    figsize = (ncolumns*4,nrows*4)
 
                plt.rcParams['axes.labelsize'] = 24
                
                howto = ["plt.rcParams['axes.labelsize'] = 24"]
                howto.append("attr = ['{:}']".format(attr))
                howto.append("xdat = x[attr]")
                howto.append("itimes = xdat.groupby('ec{:}').groups[1]".format(xdat.eventCode))
                howto.append("df = xdat.isel_points(time=itimes)")
                self.results[alias]['text'].update({'begin':{'howto': howto}})
 
                for iaxis in [1,2]:
                    df.mean(axis=iaxis).plot(robust=True)
                    
                    howto = ["df.mean(axis={:}).plot(robust=True)".format(iaxis)]
                    howto.append("plt.show()")
                    name = df.dims[iaxis].replace(alias+'_','')
                    plt_type = 'events_vs_{:}'.format(name) 
                    self.add_plot(alias, plt_type, howto=howto)
                           
                if nevents <= nmax_events:
                    #plt.rcParams['axes.labelsize'] = labelsize
                    df.plot(col='points', col_wrap=ncolumns, robust=True)
                    
                    #howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
                    howto = ["df.plot(col='points', col_wrap={:}, robust=True)".format(ncolumns)]
                    howto.append("plt.show()")
                    name = attr.replace(alias+'_','')
                    plt_type = 'events_{:}_2D'.format(name) 
                    self.add_plot(alias, plt_type, howto=howto)

# Methods to build html file

    def to_html(self, path=None, **kwargs):
        """Write out html file
        """
        self._init_html(path=path)
        self._add_html()
        self._close_html()


    def _init_html(self, path=None, **kwargs):

        self.html = output_html.report(self.exp, self.run, 
                 title=self.title,
                 css=('css/bootstrap.min.css','jumbotron-narrow.css','css/mine.css'),
                 script=('js/ie-emulation-modes-warning.js','js/jquery.min.js','js/toggler.js','js/sticky.js'),
                 output_dir=self.output_dir)

        event_times = pd.to_datetime(self._xdat.time.values)
        run_time = (event_times.max()-event_times.min()).seconds
        minutes,fracseconds = divmod(run_time,60)

        self.html.start_block('Data Summary', id="metadata")
        
        self.html.start_subblock('Data Information',id='datatime')
        sformat = 'Start Time: {:}<br/>End Time: {:}<br/>Duration: {:} seconds ({:02.0f}:{:02.0f})'
        self.html.page.p('Total events: {:}'.format(len(event_times) ) )
        self.html.page.p(sformat.format(event_times.min(), event_times.max(), run_time, minutes,fracseconds) )
        self.html.page.p('Total files: {:}<br/>Total bytes: {:} ({:0.1f} GB)<br/>'.format(
                self._nfile,self._nbytes,self._nbytes/(1000.**3)))
        self.html.page.p('Report time: {:}'.format(time.ctime()))
        self.html.end_subblock()

        self._make_PyDataSource_html()

        self.html.end_block()         
 
    def _make_PyDataSource_html(self, **kwargs):
        import PyDataSource

        self.html.start_subblock('Access the Data')
        self.html.page.p('Access event data with PyDataSource python module on a psana node:')
        #ana_env = '.  /reg/g/psdm/etc/ana_env.sh'
        idatasource = '~koglin/bin/idatasource --exp {:} --run {:}'.format(self.exp, self.run) 
        #self.html.add_textblock('\n'.join([ana_env, idatasource]))
        self.html.add_textblock('\n'.join([idatasource]))
        self.html.end_subblock()

        try:
            pyds = PyDataSource.PyDataSource
        except:
            pyds = PyDataSource

        self.html.add_textblock(pyds.__doc__, 
                subblock='HowTo access event data with PyDataSource python package', hidden=True)
        
        self.html.page.p('Analyze run summary data on a psana node using pylab, pandas and xarray:')
        ixarray = '. ~koglin/bin/ixarray --exp {:} --run {:}'.format(self.exp, self.run)
        self.html.add_textblock('\n'.join([ixarray]))
        self.html.add_textblock(str(self._xdat), 
                subblock='View of RunSummary data with xarray python package ', hidden=True)
        self.html.page.p('See http://xarray.pydata.org for details on how to use data using xarray.')
        self.html.page.p('HDF5 summary file (using netCDF format) located at:')
        self.html.add_textblock('{:}'.format(self.nc_file))
        self.html.page.p('Created by Jason Koglin:')
        self.html.page.p('For questions and feedback contact koglin@slac.stanford.edu')
 
    def _add_html(self, table_caption='Hover to highlight.', show_attrs=['attrs'], **kwargs):
        for catagory, item in self.results.items():
            self.html.start_block('{:} Data'.format(catagory), id='{:}_data'.format(catagory))
            ahowto = []
            
            if item['figure']:
                ahowto.append('# For interactive plotting -- plt.ioff() to turn off interactive plotting.')
                ahowto.append('plt.ion()')
                ahowto.append('# Alternatively make plt.show() after each plot and close window to make next')

            datatyp = 'text'
            for name, data in item[datatyp].items():
                howto = data.get('howto')
                if howto is not None:
                    howto_step = "# Howto {:} {:}:\n".format(name, catagory)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            datatyp = 'table'
            for name, data in item[datatyp].items():
                howto = data.get('howto')
                if howto is not None:
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                else:
                    howto_step = ''
                
                if data.get('format'):
                    formatters = {a: b.format for a,b in data.get('format').items()}
                    formatterstr = {a: b+'.format' for a,b in data.get('format').items()}
                else:
                    formatters = None

                df = data.get('DataFrame')
                doc = data.get('doc')
                dfname = data.get('name')
                if df is not None:
                    if name in show_attrs:
                        hidden = False
                    else:
                        hidden = True
                    
                    if formatters:
                        dfstr = df.to_string(justify='left',formatters=formatters)
                        #if dfname:
                        #    howto_step += '\n# to reprsent with formatting'
                        #    howto_step += "\nprint {:}.to_string(justify='left',formatters={:})".format(dfname, formatterstr)
                    else:
                        dfstr = str(df)
                        #if dfname:
                        #    howto_step += "\n print {:}".format(dfname)
                    
                    if howto_step:
                        ahowto.append(howto_step)

                    self.html.add_textblock(dfstr, doc=doc, 
                            subblock='{:} {:} {:}'.format(catagory,name,datatyp), 
                            hidden=hidden)

            datatyp = 'figure'
            for name, data in item[datatyp].items():
                png = data.get('png')
                if png:
                    self.html.start_subblock('{:} {:} {:}'.format(catagory, name, datatyp), hidden=True)
                    self.html.page.a( markup.oneliner.img(src=png,style='width:100%;'), 
                            href=png )
                    self.html.end_subblock(hidden=True)
            
                howto = data.get('howto')
                if howto:
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            if ahowto:
                self.html.add_textblock(ahowto, 
                        subblock='HowTo make {:} tables and figures.'.format(catagory), 
                        id='howto_{:}'.format(catagory.replace(' ','_')), 
                        hidden=True)

            self.html.end_block()         

    def _close_html(self, **kwargs):

        # this closes the left column
        self.html.page.div.close()
        
        self.html.mk_nav()
        
        self.html._finish_page()
        self.html.myprint(tofile=True)


#                Axes = pd.tools.plotting.scatter_matrix(dfscat, figsize=figsize)
#                #y ticklabels
#                [plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
#                #x ticklabels
#                [plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
#                #y labels
#                [plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]
#                #x labels
#                [plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]



