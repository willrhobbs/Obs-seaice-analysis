load "~/scripts/NCL/ncl_funcs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;read Antarctic total SIA, 
;plot climatology, with specified years

begin

  clst  = 1979 ;1989 ;climatology end
  clen  = 2008  ;climatology start
  years = (/2022,2023,2024/)
  tres  = "monthly" ; time resolution (monthly or daily)

  isDiff = False  ;if True plot dSIA/dt, plot SIA otherwise 
  MAq    = 1;21     ;smoothing period (time step, days or months)


  var  = "SIA"
  src = "cdr"
  diri = "~/data/NSIDC/"
  fili  = diri+"NSIDC_SH_totalSIA_"+tres+"_v0"+(/4,2/)+"r00.nc"
  iNam  = var+"_"+src


  diro = "plots/"
  if (isDiff) then
    filo = src+"_d"+var+"dt_"+tres+"_cycle_plot"
  else
    filo = src+"_"+var+"_"+tres+"_cycle_plot"
  end if

  oType = "png"

;*******************
;read data
;********************


  fi   = addfiles(fili, "r")

 
  SIA   = fi[:]->$iNam$
  if (isnan_ieee(SIA@_FillValue)) then
    SIA@_FillValue :=-999.
    replace_ieeenan(SIA, SIA@_FillValue,0);just in case...
  end if

  if (isDiff) then
    SIA = center_finite_diff(SIA, 1, False, 0)
    SIA = SIA * 1.e3
  end if
  

  SIA = runave(SIA, MAq, 0)



;deal with time axes

  ;change order so that 'year' starts in March
  if (tres.eq."daily") then
    ntim  = 366
 
    otime = yyyymmdd_time(2000,2000,"integer")%10000 ;MMDD
    pq       = array_append_record(ispan(60,365,1),ispan(0,59,1),0)
    otime    = otime(pq)
    YYYYMMDD = cd_calendar(time,-2)


    st    = ind(YYYYMMDD.eq.clst*10000+101)
    en    = max(ind(YYYYMMDD.le.clen*10000+1231))

    itmp = SIA(st:en)
    iMMDD = YYYYMMDD(st:en)%10000


    clim = new((/3,ntim/),typeof(SIA))
    do d = 0, 365
      dd       := ind(iMMDD.eq.otime(d))
   
      clim(1,d) = avg(itmp(dd))
      clim(0,d) = stddev(itmp(dd))
    end do

    clim(2,:)   = clim(1,:) + 1.96*clim(0,:)
    clim(0,:)   =  clim(1,:) - 1.96*clim(0,:)



    nanom      = dimsizes(years)
    anom       = new((/nanom,ntim/), typeof(SIA))
    anom!1     = "MMDD"
    anom&MMDD  = otime

    do i = 0, nanom-1

      st   = ind(YYYYMMDD.eq.years(i)*10000+301)
      en   = ind(YYYYMMDD.eq.(years(i)+1)*10000+228)
      nday = en-st
      anom(i,:nday) = (/SIA(st:en)/)
    end do

  else  ;monthly time step

    ntim  = 12
    clim  = new((/3,ntim/), typeof(SIA))
    nanom = dimsizes(years)
    anom  = new((/nanom,ntim/), typeof(SIA))

    time = fi[:]->YYYYMM

    st    = ind(time.eq.clst*100+1)
    en    = max(ind(time.le.clen*100+12))



    do m = 0, ntim-1

      clim(1,m) = avg(SIA(st+m:en:12))
      clim(0,m) = stddev(SIA(st+m:en:12))
    end do

    clim(2,:)   = clim(1,:) + 2.*clim(0,:)
    clim(0,:)   = clim(1,:) - 2.*clim(0,:)

    YYYY = time/100
   
   
    do i = 0, nanom-1

      
      ii       := ind(YYYY.eq.years(i))
      anom(i,:dimsizes(ii)-1) = (/SIA(ii)/)
    end do

    ;change order so that year starts in March
 ;   pq = array_append_record(ispan(2,11,1),(/0,1/),0)
pq = ispan(0,11,1)
    clim = clim(:,pq)
    anom = anom(:,pq)


  end if
 



  
;*********************
;plot climatology
;*********************

  wks = gsn_open_wks(oType, diro+filo)

  res                      = True
  res@gsnMaximize          = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@vpHeightF            = 0.7
  res@vpWidthF             = 1.
  
  if (isDiff) then
    res@tiYAxisString      = "dSIA/dt (10~S~6~N~ km~S~2~N~d~S~-1~N~)"
    res@gsnYRefLine        = 0. 
  else
    res@tiYAxisString      = "SIA (10~S~6~N~ km~S~2~N~)"
  end if
	
  res@trXMinF              = 0
  res@trXMaxF              = ntim-1
  res@tmXBMode             = "Explicit"
  if (tres.eq."daily") then
    res@tmXBValues         = ispan(15,366,60)
    res@tmXBLabels         = monref(otime(res@tmXBValues)/100-1)
    res@tmXBMinorValues    = ispan(45,366,60)
    res@gsnLeftString      = "("+MAq+" day mean)"
  else
    res@tmXBValues         = ispan(0,11,1)
    res@tmXBLabels         = monref(pq)
 ;   res@gsnLeftString      = "("+MAq+" month mean)"
  end if  


  res@xyLineThicknesses    = (/1.,4.,1./)
  res@xyMonoDashPattern    = True
  res@gsnXYAboveFillColors = "gray"
  res@gsnXYBelowFillColors = res@gsnXYAboveFillColors
  res@xyLineColors         = (/  res@gsnXYAboveFillColors,"black",res@gsnXYAboveFillColors/)

  climplot = gsn_csm_xy(wks, ispan(0,ntim-1,1), clim, res)

  delete([/  res@gsnXYAboveFillColors,  res@gsnXYBelowFillColors,  res@xyLineThicknesses/])

  ;specif year anomalies
  res@xyLineThicknessF     = 4.	
  res@xyLineColors        := (/"orange","red1","red4"/)
  res@pmLegendDisplayMode  = "Always"
  res@xyExplicitLegendLabels  = years
  res@pmLegendWidthF          = 0.2

  if (isDiff) then
 
    res@pmLegendParallelPosF   = 0.1
   res@pmLegendOrthogonalPosF = -0.5
  else
   res@pmLegendOrthogonalPosF = -0.5
  end if

  oplot    = gsn_csm_xy(wks, ispan(0,ntim-1,1), anom, res)
  overlay(climplot,oplot)

  draw(climplot)
  frame(wks)
  


end
