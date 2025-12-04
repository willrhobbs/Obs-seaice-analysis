#!/usr/bin/env python
# coding: utf-8

# ### Generic function to read ocean T/S data

# In[1]:


import xarray as xr
import cf_xarray as cf
import numpy as np
import glob
import os


# In[82]:


def ocean_read(src: str, expt: 'obs', var: str, start_year: int, end_year: int, 
               latmin = -90, latmax = 90, zmin = 0., zmax = 6000., freq = '1mon', chunking = -1):

    yrstring = [str(yr) for yr in  np.arange(start_year,end_year+1) ]        #string list of years
    
    #get filepath 
    if (src == 'ACCESS-OM2'):

        import intake

        catalog = intake.cat.access_nri.search(model = src,
                                              variable = var,
                                               frequency = freq
                                              )
        #get filepaths from intake catalog
        pattern = catalog[expt].search(variable = var).df['path']     

        if freq == 'fx':    #time-fixed variables
            fpath = pattern[0]
        else:
            #restrict filepath to years of interest
            fpath = sorted([f for f in pattern if any(yr in f for yr in yrstring)])

    else:
        diri = '/g/data/gv90/wrh581/'+src+'/'
        if src == 'EN4':
            pattern = os.path.join(diri, f"EN.4.2.2.?.analysis.l09.*.nc")
        elif src == 'ORAS5':
            pattern = os.path.join(diri+var, f"ORAS5_{var}_monthly_SOcean_*.nc")
        
        fpath = sorted([f for f in glob.glob(pattern) if any(yr in f for yr in yrstring)])
    


    #preprocess/spatial selection

        
    #set up preprocess
    if src == 'ORAS5':   #needs "special" treatment

        def _selection(ds):
            da = ds[var].isel(x = slice(0,1440))  # ORAS5 seems to have extra lon points/wraparounds
                    
            #vertical selection
            da = da.sel(deptht = slice(zmin, zmax))
    
            #meridional selection (noting that ORAS has 2-d lat/lon coords
            if (latmax - latmin) < 180.:
                lat = da.nav_lat.compute()
                da = da.where((lat >= latmin) & 
                              (lat >= latmax), drop = True)
    
            return da

        
    else:
 
        if freq == 'fx':  #only really the case for fx fields

            out = xr.open_dataset(fpath,decode_timedelta=False)[var]
            dims = out.dims

            if len(dims) == 3:
                    space_range = {dims[0] : slice(zmin, zmax),
                               dims[1] : slice(latmin, latmax)}
            if len(dims) == 2:
                    space_range = {dims[0] : slice(latmin, latmax)}
            
            out = out.sel(**space_range)
            
        else:
        
            def _selection(ds): 
                dims = xr.open_dataset(fpath[0], decode_timedelta=False)[var].dims
    
                if len(dims) == 4:
                    space_range = {dims[1] : slice(zmin, zmax),
                               dims[2] : slice(latmin, latmax)}
                if len(dims) == 3:
                    space_range = {dims[1] : slice(latmin, latmax)}
                
                da = ds[var].sel(**space_range)
                return da
    
            out = xr.open_mfdataset(fpath,preprocess = _selection, decode_timedelta=False, chunks = chunking, parallel=True)[var]
        
    return out

