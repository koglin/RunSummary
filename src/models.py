
def linear_model(x, attr, xvars, fit_intercept=None, name=None, 
        cut=None, residuals=True, quiet=True,
        model='LinearRegression'):
    """Make a linear model for attr based on xvars as free parameters.
       Currently only model='LinearRegression' implmented.
       Uses scikit-learn.
       residuals:  Name of attribute for residuals (default: attr+"_residuals")
    """
    if model is not 'LinearRegression':
        raise Exception("Currently only model='LinearRegression' implmented.")

    from sklearn.linear_model import LinearRegression
    import pandas as pd
    import numpy as np
    import xarray as xr
    lm = LinearRegression()  
    if not name:
        name = '{:}_model'.format(attr)

    if not quiet:
        print '\nUsing scikit-learn LinearRegression to build model for {:} from variables:\n  {:}'.format(attr, str(xvars))

    allattrs = xvars + [attr]
    if cut:
        allattrs += [cut]

    xx = x.reset_coords()[allattrs].where(np.isfinite(x.reset_coords()[attr]), drop=True)
    df_xvars0 = xx[xvars].to_dataframe()
    if cut:
        df_xvars = xx[xvars].where(xx[cut] == 1, drop=True).to_dataframe()
        xdata = xx[attr].where(xx[cut] == 1, drop=True).data
        if not quiet:
            print '\nUsing following cut in buildiing model'
            print xx[cut]
    else:
        df_xvars = df_xvars0
        xdata = xx[attr].data
    
    if fit_intercept is not None:
        lm.fit_intercept = fit_intercept
    
    lm.fit(df_xvars, xdata) 
    #ft = pd.DataFrame(zip(df_xvars.columns,lm.coef_), columns=['params','estimatedCoefficients'])
   
    x[name] = xr.DataArray(lm.predict(df_xvars0), coords=[('time', df_xvars0.index)])
    #x[name] = (['time'], lm.predict(df_xvars0))
    x[name].attrs['unit'] = x[attr].attrs.get('unit','')
    x[name].attrs['doc'] = 'LinearRegression scikit-learn model for {:} training data'.format(attr)
    x[name].attrs['model'] = model
    x[name].attrs['params'] = lm.get_params()
    x[name].attrs['variables'] = xvars
    x[name].attrs['coef_'] = lm.coef_
    x[name].attrs['intercept_'] = lm.intercept_
    x[name].attrs['score'] = lm.score(df_xvars, xdata)
    if not quiet:
        print '\n****Model Results****'
        print x.reset_coords()[name]

    if residuals:
        if not isinstance(residuals, str):
            residuals = '{:}_residuals'.format(attr)

        x[residuals] = (['time'], x[attr]-x[name])
        x[residuals].attrs['doc'] = 'Residuals for {:} based on LinearRegression model {:}'.format(attr, name)
        if not quiet:
            print '\n****Model Residuals****'
            print x.reset_coords()[residuals]

    return x


