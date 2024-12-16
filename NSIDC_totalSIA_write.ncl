
;write total SIA and SIE from NSIDC SIC data

begin

<<<<<<< HEAD
  tstep   = "monthly"
;year =     ; only used for daily data
;  version = "v04r00"
=======
  tstep   = "daily"
year = 2024    ; only used for daily data
  version = "v02r00"
>>>>>>> 029920c3801864913f2952b8a9cf2b90d74acb62
  hemi    = "SH"
;  diri    = "~/data/NSIDC/raw/"+hemi+"/"+tstep+"/"
diri = "/g/data/gv90/P6_data/NSIDC/G10016_V2/south/daily/"
;diri = "/g/data/gv90/P6_data/NSIDC/G02202_V4/south/aggregate/"

  src     = (/"cdr","nsidc_bt"/)

  ;grid cell area file (comment out if assumuning 25km equal area
  Afil    = "/g/data/jk72/wrh581/NSIDC/NSIDC0771_CellArea_PS_S25km_v1.0.nc";"~/data/NSIDC/NSIDC0771_CellArea_PS_S25km_v1.0.nc";
  aNam    = "cell_area"


  diro    = "/g/data/jk72/wrh581/NSIDC/SIE_daily/";"~/data/NSIDC/";
  filo    = "NSIDC_"+hemi+"_totalSIA_"+tstep+"_"+version+".nc"


;******************
;read data
;*****************

;for daily data, only process the listed year
;for monthly data do all files at once

  if (tstep.eq."daily") then
     fili = systemfunc("ls "+diri+"seaice_conc_"+tstep+"_*sh_"+year+"*_"+version+".nc") 	
     filo = "NSIDC_"+hemi+"_totalSIA_"+tstep+"_"+year+"_"+version+".nc"
  else
     fili = systemfunc("ls "+diri+"seaice_conc_"+tstep+"_*sh_*_"+version+".nc")
  end if

  fi   = addfiles(fili,"r")


;grid cell area, use a data file or fixed 25 km^2 value
  if (isvar("Afil")) then
    fa = addfile(Afil, "r")
    A  = fa->$aNam$
    A := tofloat(A) * 1.e-12  ;convert from m^2 to million km^2
  else
    A  = 25.^2 * 1.e-6 ;grid cell area (million km2)
  end if


;read sic for each data source
;calc SIA and and SIE

  aList = NewList("fifo")  ;list for SIA
  eList = NewList("fifo")  ;list for SIE

  nsrc  = dimsizes(src)
  do i = 0, nsrc-1 
    iNam  = src(i)+"_seaice_conc"
    if (tstep.eq."monthly") then
      iNam = iNam+"_monthly"
    end if
    sic  := fi[:]->$iNam$
    sic  := byte2flt(sic)
    sic   = where(sic.lt.0..or.sic.gt.1., sic@_FillValue, sic)

    ;calc SIA
    if (isvar("Afil")) then
      tmp  := sic * conform(sic,A,(/1,2/))
    else
      tmp  := sic * A
    end if
    ListAppend(aList,dim_sum_n(tmp,(/1,2/)))

    ;calc SIE
    sic   = where(sic.gt.0.15, 1., 0.)
    if (isvar("Afil")) then
      tmp  := sic * conform(sic,A,(/1,2/))
    else
      tmp  := sic * A
    end if
    ListAppend(eList,dim_sum_n(tmp,(/1,2/)))

  end do
  delete(sic)

;time corrections  
  time   = fi[:]->time

  if (tstep.eq."monthly") then 
    YYYYMM = cd_calendar(time,-1)
    YYYYMM@units = "YYYYMM"
    pq   = dim_pqsort(time,2)
    ;mask Dec 1987-Jan 1988
    ;should be missing but instead they're zero)
    ii     = ind(YYYYMM.eq.198712.or.YYYYMM.eq.198801)
  else
    YYYYMMDD = cd_calendar(time,-2)
    YYYYMMDD@units = "YYYYMMDD"
    pq   = dim_pqsort(time,2)
    ;mask Dec 1987-Jan 1988
    ;should be missing but instead they're zero)
    ii     = ind((YYYYMMDD/100).eq.198712.or.(YYYYMMDD/100).eq.198801)
  end if

  do i = 0, nsrc-1
    aList[i] = aList[i](pq)
    eList[i] = eList[i](pq)

    if (.not.ismissing(ii(0))) then
      aList[i](ii) = aList[i]@_FillValue
      eList[i](ii) = eList[i]@_FillValue
    end if

  end do




;****************
;write to file
;****************

  system("rm "+diro+filo)
  fo = addfile(diro+filo,"c")

  fAtt                     = True
  fAtt@description         = "Total "+hemi+" sea ice area and extent from NSIDC climate data record"
  fAtt@version             = version
  fAtt@temporal_resolution = tstep
  fAtt@time_period         = min(time)+"-"+max(time)
  fAtt@creation_date       = systemfunc("date")
  fAtt@script              = get_script_name()
  fileattdef(fo,fAtt)

  filedimdef(fo,"time",-1,True)
  filevardef(fo,"time",typeof(time),"time")
  filevarattdef(fo,"time",time)
  fo->time = (/ time /)

  if  (tstep.eq."monthly") then
    filevardef(fo,"YYYYMM","integer","time")
    filevarattdef(fo,"YYYYMM",YYYYMM)
    fo->YYYYMM = (/ YYYYMM /)
  else
    filevardef(fo,"YYYYMMDD","integer","time")
    filevarattdef(fo,"YYYYMMDD",YYYYMMDD)
    fo->YYYYMMDD = (/ YYYYMMDD /)
  end if

  vAtt            = 1.
  vAtt@long_name  = "total sea ice area"
  vAtt@units      = "10^6 km^2"
  vAtt@_FillValue = aList[0]@_FillValue

  do i = 0, nsrc-1
   
    vAtt@long_name  = str_sub_str( vAtt@long_name, "extent", "area")
    filevardef(fo,"SIA_"+src(i),"float","time")
    filevarattdef(fo,"SIA_"+src(i),vAtt)
    fo->$("SIA_"+src(i))$ = (/ aList[i] /)

    vAtt@long_name  = str_sub_str( vAtt@long_name, "area", "extent")
    filevardef(fo,"SIE_"+src(i),"float","time")
    filevarattdef(fo,"SIE_"+src(i),vAtt)
    fo->$("SIE_"+src(i))$ = (/ eList[i] /)
  end do


end
