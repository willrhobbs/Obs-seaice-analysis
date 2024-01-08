load "~/scripts/NCL/ncl_funcs.ncl"




;plot monthly timeseries of SIA for each individual month


begin

  ;climatology start/end dates
  clst  = 1979
  clen  = 2008  

  diri  = "~/data/NSIDC/"
  fili  = "NSIDC_SH_totalSIA_monthly_v0"+(/4,2/)+"r00.nc"
  iNam  = "SIA_cdr"

  diro  = "~/plots/ANTARCTIC/2021_SIE_MIN/"
  filo  = "SIA_anoms_tplot_bymonth"
  oType = "X11"



;************************
;read data
;***********************
  
  fi  = addfiles(diri+fili,"r")
  sia = fi[:]->$iNam$  


 
;split up into seperate months
;(first pad out any missing months so that each year is 'complete'

  ;start/end dates
  stmon  = sia&time(0)%100
  styear = sia&time(0)/100
  enmon  = max(sia&time)%100
  enyear = max(sia&time)/100


  ;pad first year
  if (stmon.ne.1) then  
    sia := array_append_record(new(stmon-1,typeof(sia)), sia, 0)
  end if 

  ;end month/year

  ;pad final year
  if (enmon.ne.12) then 
    sia := array_append_record(sia, new(12-enmon, typeof(sia)),0)
  end if

  sia := month_to_monyear(sia)  ;month, year
  sia&year = ispan(styear, enyear,1)


  sia = sia - conform(sia, dim_avg(sia(:,{clst:clen})), 0)


;***********************
;plot
;************************



  wks = gsn_open_wks(oType, diro+filo)
  
  
  res                         = True
  res@gsnDraw                 = True
  res@gsnFrame                = True
  res@vpHeightF               = 0.5
  res@vpWidthF                = 1.
  res@gsnMaximize             = True
  res@pmLegendDisplayMode     = "Always"

  res@gsnYRefLine             =  0.
  res@tiYAxisString           = "SIA anom (10~S~6~N~ km~S~2~N~)"
  res@xyLineThicknessF        = 2.5
  res@xyMonoDashPattern       = True
  res@xyLineColors            = read_colormap_file("srip_reanalysis") 
  


  res@xyExplicitLegendLabels  = monref
  res@lgLabelFontHeightF      = 0.01
  res@pmLegendHeightF         = 0.2
  res@pmLegendWidthF          = 0.2
  res@pmLegendOrthogonalPosF  = -0.55
  res@pmLegendParallelPosF    = 0.42
  

  res@gsnLeftString           = "("+clst+"-"+clen+" climatology)"
  plot = gsn_csm_xy(wks, sia&year, sia, res)
  

  



end