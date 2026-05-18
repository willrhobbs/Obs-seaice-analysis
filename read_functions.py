#!/usr/bin/env python
# coding: utf-8

# ### Generic function to read ocean T/S data

# In[1]:


import xarray as xr
import cf_xarray as cf
import numpy as np
import glob
import os
import intake
import pandas as pd



def ocean_read(src: str, var: str,  start_year: int, end_year: int, 
               latmin = -90, latmax = 90, zmin = 0., zmax = 6000., expt = 'obs', freq = '1mon',chunking = -1):

    yrstring = [str(yr) for yr in  np.arange(start_year,end_year+1) ]        #string list of years
    
    #get filepath 
    if (src == 'ACCESS-OM2'):

        import intake
        catalog = intake.cat.access_nri.search(model = src,
                                              variable = var,
                                              frequency = freq)

        pattern = catalog[expt].search(variable = var).df['path']     #get filepaths from intake catalog
        if freq == 'fx':
            fpath = pattern[0]
        else:
            fpath = sorted([f for f in pattern if any(yr in f for yr in yrstring)])  #restric filepath to years of interest

    else:
        diri = '/g/data/gv90/wrh581/'+src+'/'
        if src == 'EN4':
            pattern = os.path.join(diri, f"EN.4.2.2.?.analysis.l09.*.nc")
        elif src == 'IAP':
            pattern = os.path.join(diri+var, f"{src}*_{var.capitalize()}_monthly_*.nc")
        else:
            pattern = os.path.join(diri+var, f"{src}*_{var}_monthly_*.nc")
            
        fpath = sorted([f for f in glob.glob(pattern) if any(yr in f for yr in yrstring)])
    

    #preprocess
    if src == 'ORAS5':   #needs "special" treatment

        def _preprocess(ds):
            da = ds[var].isel(x = slice(0,1440))  # ORAS5 seems to have extra lon points/wraparounds
            return da 
            
    elif src == 'IAP':
                      
        def _preprocess(ds: xr.Dataset) -> xr.Dataset: 
            dt_time = pd.to_datetime(ds['time'].values, format = "%Y%m") 
            ds = ds.assign_coords(time = dt_time)
            return ds

    else:
        _preprocess = None

    out = xr.open_mfdataset(fpath,preprocess = _preprocess, chunks = chunking, parallel=False, data_vars = 'minimal', decode_timedelta=False)[var]
    if src == 'IAP': 
        out = out.chunk(chunking)

        
    #spatial selection


   # Depth selection (lazy) ----
    for zname in ("depth", "deptht", "st_ocean", "depth_std"):
        if zname in out.coords or zname in out.dims:
            out = out.sel({zname: slice(zmin, zmax)})
            break

    #latitude selection
    
    # ---- Latitude selection: 1D slice or 2D mask 


    # Fallback to known coordinate names
    def find_latname(da):
        for name in ["lat", "latitude", "yt_ocean", "nav_lat"]:
            if name in da.coords:
                return name

    latname = find_latname(out)
    if latname and (latmax - latmin) < 180.0:
        latc = out[latname]
        if latc.ndim == 1:
            out = out.sel({latname: slice(latmin, latmax)})
        else:
            latmask = (latc >= latmin) & (latc <= latmax) 
            
            # Reduce mask to 1-D keepers (still lazy)
            y_keep = latmask.any(dim='x').compute()
            x_keep = latmask.any(dim='y').compute()
            
            # Subset using indexers (small, in-memory boolean arrays)
            out = out.isel(y=y_keep, x=x_keep)

    
    return out

##################################################

###currently only works for OM2

def ice_read(src: str, expt: 'obs', var, start_year: int, end_year: int, latmax = -45., chunking = -1):

    yrstring = [str(yr) for yr in  np.arange(start_year,end_year+1) ]        #string list of years
    
    import intake
    catalog = intake.cat.access_nri.search(model = src,
                                          variable = var,
                                          frequency = '1mon')
    
    pattern = catalog[expt].search(variable = var).df['path']     #get filepaths from intake catalog
    fpath = sorted([f for f in pattern if any(yr in f for yr in yrstring)])  #restric filepath to years of interest
    
    xr.set_options(use_new_combine_kwarg_defaults=True)
    data = xr.open_mfdataset(fpath, data_vars = [var], parallel = True)[var]
    
    #correct time array
    import datetime as dt
    data['time'] = data.time.to_pandas() - dt.timedelta(hours=12)
    
    #add 1-d coords and select    
    catalog = intake.cat.access_nri.search(model = 'ACCESS-OM2', variable ='area_t')
    path = catalog[expt].search(variable = 'area_t').df['path'][0] 
    A = xr.open_dataset(path)['area_t']
    
    iNams = data.dims[1:] ; oNams = A.dims
    
    for i in [0,1]:
        data.coords[iNams[i]] = A[oNams[i]].values
        
    data = data.rename(({iNams[0]:oNams[0], 
                         iNams[1]:oNams[1]}))


    return data.sel(**{oNams[0] :slice(-90, latmax)}) 