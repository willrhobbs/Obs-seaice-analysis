
;write total SIA from NSIDC SIC data

begin

  tstep   = "monthly"
  version = "v02r00"
  hemi    = "SH"
  diri    = "~/data/NSIDC/raw/"+hemi+"/"+tstep+"/"

  src     = (/"cdr","nsidc_bt"/)

  ;grid cell area file (comment out if assumuning 25km equal area
  Afil    = "~/data/NSIDC/NSIDC0771_CellArea_PS_S25km_v1.0.nc"
  aNam    = "cell_area"


  diro    = "~/data/NSIDC/"
  filo    = "NSIDC_"+hemi+"_totalSIA_"+tstep+"_"+version+".nc"


;******************
;read data
;*****************

  fili = systemfunc("ls "+diri+"seaice_conc_monthly_*sh_*_"+version+".nc")
  fi   = addfiles(fili,"r")

;get SIA for each data source
  if (isvar("Afil")) then
    fa = addfile(Afil, "r")
    A  = fa->$aNam$
    A := tofloat(A) * 1.e-12  ;convert from m^2 to million km^2
  else
    A    = 25.^2 * 1.e-6 ;grid cell area (million km2)
  end if

  aList = NewList("fifo")
  nsrc  = dimsizes(src)
  do i = 0, nsrc-1
    
    iNam  = src(i)+"_seaice_conc_monthly"
    sic  := fi[:]->$iNam$
    sic  := byte2flt(sic)
    sic   = where(sic.lt.0..or.sic.gt.1., sic@_FillValue, sic)
    if (isvar("Afil")) then
      sic   = sic * conform(sic,A,(/1,2/))
    else
      sic   = sic * A
    end if
    ListAppend(aList,dim_sum_n(sic,(/1,2/)))

  end do
  delete(sic)

;time corrections  
  if (tstep.eq."monthly") then 
    time = cd_calendar(fi[:]->time,-1)
    time@units = "YYYYMM"
    pq   = dim_pqsort(time,2)
    ;mask Dec 1987-Jan 1988
    ;should be missing but instead they're zero)
    ii     = ind(time.eq.198712.or.time.eq.198801)
  else
    time = cd_calendar(fi[:]->time,-2)
    time@units = "YYYYMMDD"
    pq   = dim_pqsort(time,2)
    ;mask Dec 1987-Jan 1988
    ;should be missing but instead they're zero)
    ii     = ind((time/100).eq.198712.or.(time/100).eq.198801)
  end if

  do i = 0, nsrc-1
    aList[i] = aList[i](pq)
    
    if (.not.ismissing(ii(0))) then
      aList[i](ii) = aList[i]@_FillValue
    end if

  end do




;****************
;write to file
;****************

  system("rm "+diro+filo)
  fo = addfile(diro+filo,"c")

  fAtt                     = True
  fAtt@description         = "Total "+hemi+" sea ice area from NSIDC climate data record"
  fAtt@version             = version
  fAtt@temporal_resolution = tstep
  fAtt@time_period         = min(time)+"-"+max(time)
  fAtt@creation_date       = systemfunc("date")
  fAtt@script              = get_script_name()
  fileattdef(fo,fAtt)

  filedimdef(fo,"time",-1,True)
  filevardef(fo,"time","integer","time")
  filevarattdef(fo,"time",time)
  fo->time = (/ time /)


  vAtt            = 1.
  vAtt@long_name  = "total sea ice area"
  vAtt@units      = "10^6 km^2"
  vAtt@_FillValue = aList[0]@_FillValue

  do i = 0, nsrc-1
   
    filevardef(fo,src(i),"float","time")
    filevarattdef(fo,src(i),vAtt)
    fo->$src(i)$ = (/ aList[i] /)

  end do


end
