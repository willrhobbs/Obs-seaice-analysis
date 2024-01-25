load "~/scripts/NCL/ncl_funcs.ncl"
load "~/scripts/NCL/plot_funcs.ncl"

;plot areal-mean ocean anomalies from RG gridded ARGO data




begin

  ;spatial domain
  latmax = -55.
  latmin = -65.
  zmin   = 0.
  zmax   = 600.

  
  ;time domain
  stdat  = 200401
  endat  = 202312
  clst   = 2004
  clen   = 2018

  refdats  = (/201608/)  ;highlight specific dates (YYYYMM format)


  ;input file 
  diri   = "~/data/RG_ARGO/"
;  var    = "Salinity";"Temperature";
  vNam    = "ARGO_"+str_upper(var)+"_ANOMALY"
  fili   = "RG_ArgoClim_"+var+"_200401-"+endat+".nc"

  ;output file
  pType  = "png"
  dirp   = "plots/"
  filp   = "ARGO_"+var+"_depthtime_SOcean_"+zmin+"-"+zmax+"m"


;*******************
;read data
;*******************


  fi    = addfile(diri+fili,"r")
  time  = cd_calendar(fi->TIME,-1)



  tmp  = fi->$vNam$(:,{zmin:zmax},{latmin:latmax},:)


  dlat = 1.
  dy   = yresCalc(dlat)
  dlon = 1.
  dx   = xresCalc(dlon, tmp&LATITUDE)
  A    = conform(tmp(0,0,:,:), dx*dy, 0)
  
  
  varmn = tmp(:,:,0,0) 
  varmn = wgt_areaave2(tmp, A, 0)  ;get areal mean

  


  ;get anomaly 
  st    = ind(time.eq.clst*100+1)
  en    = ind(time.eq.clen*100+12)

  varmn&TIME := time
  anom      = varmn({stdat:endat},:)
  sig       = anom                ;array for st. deviation
  nmon      = min((/12,dimsizes(anom(:,0))/))

  


  do m = 0, nmon-1
    mon           = stdat%100 + m        
    clim          = dim_avg_n(varmn({clst*100+mon:clen*100+mon:12},:),0)
    anom(m::12,:) = anom(m::12,:) - conform( anom(m::12,:),clim,1)

    clim          = dim_stddev_n(varmn({clst*100+mon:clen*100+mon:12},:),0)
    sig(m::12,:)  =  conform(sig(m::12,:),clim,1)
  end do
    
  delete([/clim,tmp,varmn, A/])
  

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
  res@trYReverse           = True
  

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
    res@cnMinLevelValF       = -0.5
    res@cnFillPalette        = "cmocean_balance"
  else
    res@cnMinLevelValF       = -0.05
    res@cnFillPalette        = "cmocean_delta"
  end if

  res@cnMaxLevelValF       = -res@cnMinLevelValF
  res@cnLevelSpacingF      = res@cnMaxLevelValF/10.
  res@lbLabelBarOn         = True
  res@cnInfoLabelOn        = False
  res@cnLineLabelsOn       = False



  res@gsnLeftString       = (-latmin)+"-"+(-latmax)+"~S~o~N~S"
  res@gsnRightString       = "("+clst+"-"+clen+" clim)"
  res@tiMainString         = var




  time     := anom&TIME
  anom&TIME = tdum

  ;mask out near-zero anoms
  anom      = where(abs(anom).lt.res@cnLevelSpacingF,anom@_FillValue, anom )

  plot =  gsn_csm_contour(wks,anom(PRESSURE|:,TIME|:), res)


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

  oplot = gsn_csm_contour(wks, sig(PRESSURE|:,TIME|:), res)

  overlay(plot,oplot)

;add reference dates
  nref                  = dimsizes(refdats)
  lines                 = new((/nref/), graphic)
  resl                  = True
  resl@gsLineThicknessF = 5.
  resl@gsLineColor      = "red"

  do i = 0, nref-1
    xref = ind(time.eq.refdats(i))  
    lines(i) = gsn_add_polyline(wks,plot,(/xref,xref/),(/0.,zmax/), resl)
  end do




  draw(plot)
  frame(wks)


end
