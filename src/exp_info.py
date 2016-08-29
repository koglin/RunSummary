
from RegDB import experiment_info
from glob import glob
import re
import operator
import sys
import os
import time
import traceback
from LogBook.runtables import RunTables


class ExperimentInfo(object):

    def __init__(self, exp=None, instrument=None, station=0, exper_id=None,
                exp_dir=None, xtc_dir=None, h5_dir=None):
        if exp:
            self.exp = exp
            self.exper_id = experiment_info.name2id(exp)
        
        elif exper_id:
            self.exper_id = exper_id
            self.exp = experiment_info.getexp(exper_id)
        
        elif instrument:
            self.exp = experiment_info.active_experiment(instrument, station=station)
            self.exper_id = experiment_info.name2id(exp)

        if not instrument:
            instrument = self.exp[0:3]
            
        self.station = station
        self.instrument = instrument

        if not exp_dir:
            exp_dir = "/reg/d/psdm/{:}/{:}".format(self.instrument, self.exp)

        if not xtc_dir:
            xtc_dir =  "{:}/xtc".format(exp_dir, self.exp)

        if not xtc_dir:
            h5_dir =  "{:}/hdf5".format(self.exp_dir, self.exp)
        
        self.exp_dir = exp_dir
        self.h5_dir = h5_dir
        self.xtc_dir = xtc_dir

        # Setup elog pswwww RunTables
        self._RunTables = RunTables(**{'web-service-url': 'https://pswww.slac.stanford.edu/ws-kerb'})

        self._user_tables = self._RunTables.usertables(exper_name=self.exp)
        #self.add_user_run_table(RunSummary='Run Summary')

#    def add_user_run_table(self, **kwargs):
#        """Add a user User RunTable from pswww elog server.
#        """
#        for alias, name in kwargs.items():
#            tbl = self._RunTables.findUserTable(exper_name=self.exp, table_name=alias)
#            self._user_tables.update({alias: tbl}) 
#            setattr(self, alias, tbl)

    def detectors(self, run):
        """Return a list of detector names configured in the DAQ system for the input run number.
        """
        return experiment_info.detectors(self.instrument, self.exp, run)

    @property
    def calibration_runs(self):
        return experiment_info.calibration_runs(self.instrument, self.exp)

    @property
    def runs(self):
        """Experiment run information from MySQL database and xtc directory.
        """
        if experiment_info.name2id(self.exp):
            runs_list =  experiment_info.experiment_runs(self.instrument.upper(),self.exp)
            for item in runs_list:
                runnum = item['num']
                item['xtc_files'] = glob('{:}/*-r{:04d}*.xtc'.format(
                                        self.xtc_dir,runnum))
                item['h5_files'] = glob('{:}/*-r{:04d}*.h5'.format(
                                        self.h5_dir,runnum))
        else:
            runs_list = []

        return runs_list

    @property
    def open_files(self, run=None):
        """Return a list of files created (by the DAQ system).  
           Current run if no run is specified.
        """
        return experiment_info.get_open_files(self.exper_id,run)

    def load_run_summary(self):
        """Load MySQL database experiment run summary information into a dictionary.
        """
        vrun_attrs = {}
        print 'Loading summary of {:} runs for {:} from SQL database'.format( \
                len(self.runs),self.exp)
        print 'Estimate loading time ~{:} sec'.format(len(self.runs)/4)
        for run in range(1,self.runs[-1]['num']+1):
            run_attr = experiment_info.run_attributes(self.instrument,self.exp,run)
            for a in run_attr:
                if a['name'] not in vrun_attrs:
                    vrun_attrs[a['name']] = {'class': a['class'], 'desc': a['descr'],
                                             'type': a['type'], 'val':
                                             [None for i in range(1,run)]}
                vrun_attrs[a['name']]['val'].append(a['val'])
        self.run_summary = vrun_attrs


    def show_runs(self,start=0,end=99999999,csv=False):
        """Show run summary for current experiment.
        """
        if csv:
            print '{:>7}, {:>10}, {:>8}, {:>10}, {:3}, {:2}'.format('Run',
                                'Day', 'Time', 'Length', 'xtc', 'h5')

        else:
            print '='*72
            print 'Experiment {:}'.format(self.exp)
            print '  xtc dir {:}'.format(self.xtc_dir)
            print '  hdf5 dir {:}'.format(self.h5_dir)
            print '-'*72
            print '{:>7} {:>10} {:>8} {:>10} {:3} {:2}'.format('Run', 'Day', 'Time',
                                                  'Length', 'xtc', 'h5')
            print '-'*72

        for item in self.runs:
            run = item['num']
            if run >= start and run <= end:
                datestr = time.strftime('%Y-%m-%d',
                                        time.localtime(item['begin_time_unix']))
                timestr = time.strftime('%H:%M:%S',
                                        time.localtime(item['begin_time_unix']))
                if len(item['xtc_files']) > 0:
                    xtc = 'xtc'
                else:
                    xtc = ''

                if len(item['h5_files']) > 0:
                    h5 = 'h5'
                else:
                    h5 = ''

                begin_time = item['begin_time_unix']
                end_time = item['end_time_unix']
                if end_time:
                    dtime = end_time - begin_time
                    flag = ' '
                else:
                    dtime = time.time() - begin_time
                    flag = '*'

                dmin = int(dtime/60)
                dsec = int(dtime % 60)
                if dmin > 0:
                    dtstr = '{:4}m {:02}s'.format(dmin,dsec)
                else:
                    dtstr = '{:02}s'.format(dsec)

                if csv:
                    print '{:7}, {:10}, {:8}, {:>10}, {:3}, {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5)
                else:
                    print '{:7} {:10} {:8} {:>10} {:3} {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5)

                if flag in '*':
                    print '* Currently Acquiring Data for Run {:}'.format(run)


