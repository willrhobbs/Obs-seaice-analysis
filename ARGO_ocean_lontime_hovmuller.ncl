load "~/scripts/NCL/ncl_funcs.ncl"
load "~/scripts/NCL/plot_funcs.ncl"

;plot depth-mean ocean anomalies from RG gridded ARGO data




begin

  ;spatial domain
  latmax = -50.
  latmin = -60.
  zmin   = 100.
  zmax   = 150.


 
  
  ;time domain
  stdat  = 200401
  endat  = 202404
  clst   = 2004
  clen   = 2018

  refdats  = (/201112,201608/)  ;highlight specific dates (YYYYMM format)


  ;input file 
  diri   = "~/data/RG_ARGO/"
;  var    = "Salinity";"Temperature";
  vNam    = "ARGO_"+str_upper(var)+"_ANOMALY"
  fili   = "RG_ArgoClim_"+var+"_200401-"+endat+".nc"

  ;output file
  pType  = "png"
  dirp   = "plots/"
  filp   = "ARGO_"+var+"_lontime_SOcean_"+zmin+"-"+zmax+"m"


;*******************
;read data
;*******************


  fi    = addfile(diri+fili,"r")
  time  = cd_calendar(fi->TIME,-1)



  tmp  = fi->$vNam$(:,{zmin:zmax},{latmin:latmax},:)
  tmp := dim_avg_n_Wrap(tmp, 2)  ;get merid mean

  
  
;weighted depth mean
; 

  z      = fi->PRESSURE   ;full depth range
  dz     = tmp&PRESSURE  
 
  z0     = ind(z.eq.dz(0))
  z1     = ind(z.eq.max(dz))

  nz     = dimsizes(tmp(0,:,0))
  bnds   = new((/2,nz/), typeof(z))

  bnds(1,:)  = dim_avg_n((/z(z0+1:z1+1), z(z0:z1)/),0)
  bnds(0,1:) = bnds(1,:nz-2)

  if (zmin.eq.0) then
    bnds(:,0) = (/0.,2. * z(0)/)
    bnds(0,1) = bnds(1,0)
  else
    bnds(0,0) = avg(z(z0-1:z0))
  end if

  dz = bnds(1,:) - bnds(0,:)

  tmp := dim_avg_wgt_n_Wrap(tmp, dz, 1, 1)

  delete([/dz,bnds,z/])





  ;get anomaly 
  st    = ind(time.eq.clst*100+1)
  en    = ind(time.eq.clen*100+12)

  tmp&TIME := time
  anom      = tmp({stdat:endat},:)
  sig       = anom                ;array for st. deviation
  nmon      = min((/12,dimsizes(anom(:,0))/))

  


  do m = 0, nmon-1
    mon           = stdat%100 + m        
    clim          = dim_avg_n(tmp({clst*100+mon:clen*100+mon:12},:),0)
    anom(m::12,:) = anom(m::12,:) - conform( anom(m::12,:),clim,1)

    clim          = dim_stddev_n(tmp({clst*100+mon:clen*100+mon:12},:),0)
    sig(m::12,:)  =  conform(sig(m::12,:),clim,1)
  end do
    
  delete([/clim,tmp/])
  

;correct longitude
  lon = anom&LONGITUDE
  lon = where(lon.gt.360., lon-360.,lon)
  pq  = dim_pqsort(lon, 2)
  anom = anom(:,pq)
  lon@long_name  = "longitude"
  lon@units      = "degrees_east"
  anom&LONGITUDE:= lon
  sig            = sig(:,pq)


;turn standard deviation into a signifance mask
  sig = where(abs(anom).ge.(1.96*sig), 1., 0.)

;*********************
;plot 
;**********************

  wks = gsn_open_wks(pType, dirp+filp)

  
  res                    = True
  res@gsnDraw            = False
  res@gsnFrame           = False
  res@gsnMaximize         = True
  res@vpHeightF          = 0.5
  res@vpWidthF           = 1.
  

  res@tmXBMode             = "Explicit"
  step                     = 36
  ntim                     = dimsizes(anom(:,0))
  tdum                     = ispan(1,ntim,1)
  res@tmXBValues           = tdum(::step)
  res@tmXBLabels           = anom&TIME(::step)/100
  res@tmXBMinorValues      = tdum(::12)

  res@cnLinesOn            = False
  res@cnFillOn             = True
  res@cnFillMode           = "RasterFill"
  res@cnLevelSelectionMode = "ManualLevels"

  if (var.eq."Temperature") then
    res@cnMinLevelValF       = -1.
    res@cnFillPalette        = "cmocean_balance"
  else
    res@cnMinLevelValF       = -0.1
    res@cnFillPalette        = "cmocean_delta"
  end if

  res@cnMaxLevelValF       = -res@cnMinLevelValF
  res@cnLevelSpacingF      = res@cnMaxLevelValF/10.
  res@lbLabelBarOn         = True
  res@cnInfoLabelOn        = False
  res@cnLineLabelsOn       = False


  res@gsnLeftString        = zmin+"-"+zmax+"m"
  res@gsnCenterString       = latmin+"-"+latmax+"~S~o~N~N"
  res@gsnRightString       = "("+clst+"-"+clen+" clim)"
  res@tiMainString         = var


  time     := anom&TIME
  anom&TIME = tdum

  ;mask out near-zero anoms
  anom      = where(abs(anom).lt.res@cnLevelSpacingF,anom@_FillValue, anom )

  plot =  gsn_csm_contour(wks,anom(LONGITUDE|:,TIME|:), res)


;significance overlay respources

  res@gsnLeftString        = ""
  res@gsnRightString       = ""
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = 0.5     
  res@cnLinesOn            = True
  res@cnFillOn             = False
  res@cnLineColor          = "black"
  res@cnLineThicknessF     = 4.
 ; res@cnFillColors         = (/"transparent","black"/)
 ; res@cnFillPattern        = 5


  copy_VarCoords(anom,sig)

  oplot = gsn_csm_contour(wks, sig(LONGITUDE|:,TIME|:), res)

  overlay(plot,oplot)

;add reference dates
  nref                  = dimsizes(refdats)
  lines                 = new((/nref/), graphic)
  resl                  = True
  resl@gsLineThicknessF = 4.
  resl@gsLineColor      = "black"

  do i = 0, nref-1
    xref = ind(time.eq.refdats(i))  
    lines(i) = gsn_add_polyline(wks,plot,(/xref,xref/),(/0.,360./), resl)
  end do




  draw(plot)
  frame(wks)


end
