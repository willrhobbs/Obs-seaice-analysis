#!/usr/bin/env python
# coding: utf-8

# ### Generic function to read ocean T/S data

# In[1]:


import xarray as xr
import numpy as np
import glob
import os


# In[82]:


def ocean_read(src: str, var: str, start_year: int, end_year: int, space_range: dict = {}, chunking = -1):

    yrstring = [str(yr) for yr in  np.arange(start_year,end_year+1) ]        #string list of years
    
    #get filepath 
    if (src == 'ACCESS-OM2'):

        import intake
        catalog = intake.cat.access_nri.search(model = src,
                                              variable = var,
                                              frequency = '1mon')

        pattern = catalog[expt].search(variable = var).df['path']     #get filepaths from intake catalog
        fpath = sorted([f for f in pattern if any(yr in f for yr in yrstring)])  #restric filepath to years of interest

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
            if 'vertical' in space_range:
                da = da.sel(deptht = space_range['vertical'])
    
            #horizontal selection (noting that ORAS has 2-d lat/lon coords
            if 'latitude' in space_range:
                lat = da.nav_lat.compute()
                da = da.where((lat >= space_sel['latitude'].start) & 
                              (lat >= space_sel['latitude'].stop), drop = True)
    
            if 'longitude' in space_range:
                lon = da.nav_lon.compute()
                da = da.where((lon >= space_sel['longitude'].start) & 
                               (lon >= space_sel['longitude'].stop), drop = True)
            return da

        
    else:
        import cf_xarray as cf
        def _selection(ds): 
            da = ds[var].cf.sel(**space_range)
            return da

    data = xr.open_mfdataset(fpath,preprocess = _selection, chunks = chunking, parallel=True)[var]
    
    return data

