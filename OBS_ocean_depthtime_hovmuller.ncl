load "~/scripts/NCL/ncl_funcs.ncl"
load "~/scripts/NCL/plot_funcs.ncl"

;plot areal-mean ocean anomalies from RG gridded ARGO data




begin

  ;spatial domain
  sNam   = "Weddell" ; "SOcean" 
  latmax = -55.
  latmin = -65.
  lonmin = 330.; 0.
  lonmax = 20. ;360.
  zmin   = 0.
  zmax   = 1000.

  
  ;time domain
  stdat  = 202001
  endat  = 202504
  clst   = 2006
  clen   = 2018

  refdats  = (/201609,202302/)  ;highlight specific dates (YYYYMM format)


  ;input file 
;src = "EN4"
;  diri   = "/g/data/gv90/wrh581/EN4/"   ;"~/data/RG_ARGO/"
;  var    = "salinity";"temperature";
;  vNam    = var;"ARGO_"+str_upper(var)+"_ANOMALY"
;  fili = systemfunc("ls "+diri+"EN.4.2.2.?.analysis.l09.;;{"+str_join(yyyymm_time(stdat/100,endat/100,"integer"),",")+"}.nc 2>/dev/null")
;dlat = 1.
;dlon = 1.

src = "ORAS5"
;  var    = "vosaline";"votemper";
  diri   = "/g/data/gv90/wrh581/"+src+"/"+var+"/"   
  vNam    = var
  fili = systemfunc("ls "+diri+"ORAS5_"+var+"_monthly_SOcean_{"+str_join(ispan(stdat/100,endat/100,1),",")+"}.nc 2>/dev/null")
dlon = 0.25
dlat = 0.12


 
;"RG_ArgoClim_"+var+"_200401-"+endat+".nc"

  ;output file
  pType  = "png"
  dirp   = "plots/"
  filp   = src+"_"+var+"_depthtime_"+sNam+"_"+zmin+"-"+zmax+"m"


;*******************
;read data
;*******************


  fi    = addfiles(fili,"r")
  dNams = getfilevardims(fi[0],var)

  time  = cd_calendar(fi[:]->$dNams(0)$,-1)

lon = fi[0]->nav_lon(0,:1439)
print(lon)
print(isMonotonic(lon))
exit

  if (src.eq."ORAS5") then    ; tripolar grid
     tmp  = fi[:]->$vNam$(:,{zmin:zmax},:,:1439)
    tmp!2 = "lat"
    tmp&lat := fi[0]->nav_lat(:,0)
    tmp!3 = "lon"
    tmp&lon := fi[0]->nav_lon(0,:1439)
    tmp := tmp(:,:,{latmin:latmax},:)

tmp = dim_pqsor


  else
      ;rectilinear grid
       tmp  = fi[:]->$vNam$(:,{zmin:zmax},{latmin:latmax},:)
  end if




  if (lonmin.lt.lonmax) then
    tmp := tmp(:,:,:,{lonmin:lonmax})
  else
    tmp(:,:,:,{lonmax:lonmin}) = tmp@_FillValue
  end if





  dy   = yresCalc(dlat)
  dx   = xresCalc(dlon, tmp&$tmp!2$)
  A    = conform(tmp(0,0,:,:), dx*dy, 0)
  
  
  varmn = tmp(:,:,0,0) 
  varmn = wgt_areaave2(tmp, A, 0)  ;get areal mean

  


  ;get anomaly 
  st    = ind(time.eq.clst*100+1)
  en    = ind(time.eq.clen*100+12)

  varmn&$varmn!0$ := time
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
  res@tmXBLabels           = anom&$anom!0$(::step)/100
  res@tmXBMinorValues      = tdum(::12)

  res@cnLinesOn            = False
  res@cnFillOn             = True
  res@cnFillMode           = "RasterFill"
  res@cnLevelSelectionMode = "ManualLevels"

  if (any(var.eq.(/"temperature","Temperature"/))) then
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



  res@gsnLeftString       = (-latmin)+"-"+(-latmax)+"~S~o~N~S, "+lonmin+"-"+lonmax+"~S~o~N~E"
  res@gsnRightString       = "("+clst+"-"+clen+" clim)"
  res@tiMainString         = var




  time     := anom&$anom!0$
  anom&$anom!0$ = tdum

  ;mask out near-zero anoms
  anom      = where(abs(anom).lt.res@cnLevelSpacingF,anom@_FillValue, anom )

  plot =  gsn_csm_contour(wks,anom($anom!1$|:,$anom!0$|:), res)


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

  oplot = gsn_csm_contour(wks, sig($anom!1$|:,$anom!0$|:), res)

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
