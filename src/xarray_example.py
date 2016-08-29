from pylab import *
import xarray as xr
import pandas as pd

import build_html
import psxarray

file = '/reg/d/psdm/{:}/{:}/scratch/nc/run{:04}.nc'.format(xdat.instrument,xdat.experiment,int(xdat.run[0]))
xdat.to_netcdf(file, engine='h5netcdf')

ds = PyDataSource.DataSource(exp='mfxc0116',run=30)
# Add SampleCamera detector module and make roi
roi = ((130, 160), (825, 854))
path = '/reg/neh/home/koglin/psana/current/PyDataSource/src/detectors/'
ds.add_detector('SampleCamera', module='camera', path=path, roi=roi)
# Make xarray
x = psxarray.get_xdat(ds)

file_name = '/reg/d/psdm/cxi/cxig8915/scratch/nc/run0293.nc'
x = xarray.open_dataset(file_name, engine='h5netcdf')
b = build_html.Build_html(x)
b.add_summary(variables=['CsPadKalpha_img_mean'])
attrs = ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX', 'EBeam_ebeamLTUPosX', 'EBeam_ebeamLTUAngX', 'XCS_IPM_02_xpos', 'CsPadKalpha_xrays']
b.add_detector(alias='Xcorrelation',attrs=attrs)
b.add_detector(alias='Correlate',attrs=['CsPadKalpha_xrays','FEEGasDetEnergy_f_11_ENRC','EBeam_ebeamPhotonEnergy','EBeam_ebeamCharge','EBeam_ebeamL3Energy','EBeam_ebeamEnergyBC1','EBeam_ebeamPkCurrBC1'])
b.add_html()
b.close()

#-------------
#. ~koglin/bin/ixarray --exp cxik8816 --run 28
#%time a = x.DsaCsPad_calib.sum(axis=(1,2,3))
attrs = ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX', 'FEEGasDetEnergy_f_11_ENRC', 'DsaCsPad_totalADU', 'EBeam_ebeamL3Energy']
b.add_detector(alias='Xcorrelation',attrs=attrs)
b.add_summary(variables=['DsaCsPad_image_mean','DsaCsPad_image_max'])
b.add_all()
b.add_html()
b.close()


#--------------


from RunSummary import psxarray


### MFX Wire Scan data 

from pandas.tools.plotting import scatter_matrix
import xarray

path = '/reg/d/psdm/mfx/mfxc0116/scratch/nc/'
run = 31
filename = '{:}run{:04}.nc'.format(path,run)

x = xarray.open_dataset(filename,engine='h5netcdf')

# plot 2D images for each step
x.SampleCamera_img_mean[70:100,0].plot(col='steps', col_wrap=5)

# plot 2D images for each motor position (i.e., at each step but with motor position instead of step number).
x.SampleCamera_img_mean.squeeze('codes').groupby('MFX_USR_MMN_01').mean(axis=0)[70:100].plot(col='MFX_USR_MMN_01', col_wrap=5)

# Show correlation of wire scan position and projection on y axis
x.SampleCamera_img_mean.squeeze('codes').groupby('MFX_USR_MMN_01').mean(axis=0).mean(axis=1).plot()


# Scatter Matrix
z = x[['EBeam_ebeamCharge', 'EBeam_ebeamPhotonEnergy', 'FEEGasDetEnergy_f_11_ENRC', 'FEEGasDetEnergy_f_21_ENRC']].to_array().to_pandas().T

scatter_matrix(z, alpha=0.2, figsize=(6, 6), diagonal='kde')

x.SampleCamera_img.groupby('step').mean(dim='time').mean(axis=2).mean(axis=1).plot()


#temp = 15 + 8 * np.random.randn(2, 2, 3)
#precip = 10 * np.random.rand(2, 2, 3)
#lon = [[-99.83, -99.32], [-99.79, -99.23]]
#lat = [[42.25, 42.21], [42.63, 42.59]]
#
#ds = xr.Dataset({'temperature': (['x', 'y', 'time'],  temp),
#                 'precipitation': (['x', 'y', 'time'], precip)},
#                 coords={'lon': (['x', 'y'], lon),
#                         'lat': (['x', 'y'], lat),
#                         'time': pd.date_range('2014-09-06', periods=3),
#                          'reference_time': pd.Timestamp('2014-09-05')})
#
#file = '/reg/neh/home/koglin/data/text.cn'


