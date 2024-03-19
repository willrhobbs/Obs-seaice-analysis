load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "~/scripts/NCL/ncl_funcs.ncl"


;write integrated SH ice area by longitude



begin

;zonal/mred. resolutions of output grid (degrees)
  xres  = 1.
  yres  = 0.5 


  src  = "cdr"
  ;version = "v02r00"
  pole = "SH"
  diri  = "~/data/NSIDC/raw/"+pole+"/monthly/"
  
  ;grid cell area file (comment out if assumuning 25km equal area
  Afil    = "~/data/NSIDC/NSIDC0771_CellArea_PS_S25km_v1.0.nc"
  aNam    = "cell_area"

 if (pole.eq."SH") then
 
   ymin = -87.5
   ymax = -30.5
   nlat = 58
   fili = systemfunc("ls "+diri+"seaice_conc_monthly_*sh_*_"+version+".nc")
 else
   ymin = 45.5
   ymax = 89.5
   nlat = 45
   fili = systemfunc("ls "+diri+"seaice_conc_monthly_nh_*.nc")
 end if



 ;output file name
  diro = "~/data/NSIDC/"




;***********************************
;time coords
;***********************

 print("reading time data")
 
 setfileoption("nc","SuppressClose",False) 
 
 f = addfiles(fili,"r")

 ttmp = f[:]->time

 ;sort time order
 if (isMonotonic(ttmp).lt.1) then
   pq   = dim_pqsort(ttmp,1) 
   f    = addfiles(fili(pq),"r")
   ttmp = f[:]->time
 end if   


 time = cd_calendar(ttmp,-1)

 

;******************************
;read ice data
;******************************

  print("reading ice data")


  vNam = src+"_seaice_conc_monthly"
  tmp = f[:]->$vNam$



  dimz = dimsizes(tmp)
  ntim = dimz(0)
  ny   = dimz(1)
  nx   = dimz(2)

 


  sic  = new((/ntim,ny,nx/),"float")  ;new array for data
  sic = byte2flt(tmp)

  sic   = where(sic.lt.0..or.sic.gt.1., sic@_FillValue, sic)



;***********************
;regrid
;***************************


  print("regridding")


  if (version.eq."v02r00") then
    ftmp  = systemfunc("ls "+diri+"seaice_conc_monthly_sh_197811_??????_v04r00.nc")
    ftmp := addfile(ftmp,"r")
    lat2d = ftmp->latitude
    lon2d = ftmp->longitude 
  else
    lat2d = f[0]->latitude
    lon2d = f[0]->longitude
  end if

  lat1d = ndtooned(lat2d)
  lon1d = ndtooned(lon2d)
  lon1d = where(lon1d.lt.0,360.+lon1d,lon1d)

  sic1D = reshape(sic, (/ntim,ny*nx/))

  
; Retrieve coordinates that are not missing.
; We need these later to reconstruct grid.
;
  Index = ind(.not.ismissing(lon1d))    


  gOpt                = True 
  gOpt@SrcFileName    = "NSIDC_"+yres+"_x_"+yres+"_regrid_src.nc"
  gOpt@DstFileName    = "NSIDC_"+yres+"_x_"+yres+"_regrid_dst.nc"
  gOpt@WgtFileName    = "NSIDC_"+yres+"_x_"+yres+"_regrid_wgtfile.nc"
  gOpt@ForceOverwrite = True
  ;gOpt@Debug         = True

  gOpt@SrcGridLat     = lat1d(Index)    ; Source grid
  gOpt@SrcGridLon     = lon1d(Index)
  gOptSrcGridType     = "unstructured"    
  gOpt@InterpMethod   = "bilinear"
  gOpt@SrcMask2D      = where(ismissing(sic(0,:,:)),0,1) ;land mask

;Destination grids
  nlat                = toint((ymax-ymin)/yres) + 1
  gOpt@DstGridLat     = fspan(ymin,ymax,nlat)
  gOpt@DstGridLon     = fspan(xres/2.,360.-xres/2.,toint(360/xres))

  
  if (isfilepresent(gOpt@WgtFileName)) then
    ice3D = ESMF_regrid_with_weights(sic1D(:,Index),gOpt@WgtFileName,gOpt)
  else
    ice3D = ESMF_regrid(sic1D(:,Index), gOpt)
  end if



  ice3D = where(ice3D.gt.1.,1.,ice3D)
  ice3D = where(ice3D.lt.0., 0.,ice3D)


    ;get areal cover
  dy  = yresCalc(yres)   ; merid. length (scalar: km)
  dx  = xresCalc(xres,  gOpt@DstGridLat) ;zonal length scale (nlat, km)






;integrate along meridians
  A   := dy * conform(ice3D, dx, 1)
  delete([/dx,dy/])



  ;sea ice area
  ice3D = where(ice3D.ge.0.15, ice3D, 0.)
  SIA   = dim_sum_n(ice3D * A, 1)  
  SIA@long_name = "sea ice area"
  SIA@units     = "km**2"


  ;sea ice extent
  ice3D = where(ice3D.ge.0.15, 1., 0.)
  SIE   = dim_sum_n(ice3D * A, 1) 
  SIE@long_name = "sea ice extent"
  SIE@units     = "km**2"

  ;use new 1-d lon array   
  lon = dble2flt(ice3D&lon)

  delete([/ice3D, sic,sic1D,Index/])





;*****************
;write to file
;*****************

  filo = "NSIDC_"+pole+"_icearea_"+version+"_"+time(0)+"-"+max(time)+".nc"
  system("rm "+diro+filo)
  fo = addfile(diro+filo, "c")



  fAtt               = True
  fAtt@creation_date = systemfunc("date")
  fAtt@description   = "NSIDC "+pole+" area and extent"
  fAtt@algorithm     = src
  fAtt@version       = version
  fAtt@script        = get_script_name()
  fileattdef(fo, fAtt)

  filedimdef(fo, (/"time","lon"/),(/-1,dimsizes(lon)/),(/True,False/))

  
  filevardef(fo, "time","integer","time")
  filevardef(fo, "lon",typeof(lon), "lon")
  filevarattdef(fo,"time",time)
  filevarattdef(fo, "lon",lon)


  fo->time = (/time/)
  fo->lon  = (/lon/)

  filevardef(fo,"SIA",typeof(SIA), (/"time","lon"/))
  filevarattdef(fo,"SIA",SIA)
  fo->SIA  = (/SIA/)

  filevardef(fo,"SIE",typeof(SIA), (/"time","lon"/))
  filevarattdef(fo,"SIE",SIE)
  fo->SIE  = (/SIE/)


end
