load "~/scripts/NCL/ncl_funcs.ncl"
load "~/scripts/NCL/plot_funcs.ncl"

;scatterplot of annula-mean Arctic vs Antarctic SIA


begin

  yrst = 1979
  yren = 2024
  clst = 1979  ;climatology end
  clen = 2018  ;climatology start

     ;

  diri = "~/data/NSIDC/"
  fils = "NSIDC_SH_totalSIA_monthly_"+(/"v05r00.nc","v03r00.nc"/);(/"v04r00.nc","v02r00.nc"/)
  sNam = "SIE_cdr"

  filn = "NSIDC_NH_totalSIA_monthly_v03r00.txt"

  diro = "plots/"
  filo = "NSIDC_Arctic_vs_Antarctic"
  oType = "png"



;*******************
;read SH data
;********************



  fi = addfiles(diri+fils, "r")
  YYYYMM = fi[:]->YYYYMM

  st       = ind(fi[:]->YYYYMM.eq.yrst*100+1)
  en       = ind(fi[:]->YYYYMM.eq.yren*100+12)
  tmp      = fi[:]->$sNam$(st:en)
  tmp&time := YYYYMM(st:en)

  
;annual mean
  year     = ispan(yrst, yren, 1)
  SIA_SH      = new(dimsizes(year), typeof(tmp))
  SIA_SH!0    = "year"
  SIA_SH&year = year
  SIA_SH@long_name = "Antarctic"
  
  
  do y = yrst, yren
    SIA_SH({y}) = avg(tmp({y*100+1:y*100+12}))
  end do

  SIA_SH({1987:1988}) = SIA_SH@_FillValue

;*******************
;read NH data
;********************

  SIA_NH = SIA_SH
  SIA_NH@long_name = "Arctic"

  tmp := readAsciiTable(diri+filn, 13, typeof(SIA_NH), 1)

  do y = yrst, yren
    ii         = ind(tmp(:,0).eq.y)
    SIA_NH({y}) = avg(tmp(ii,1:))
  end do

  SIA_NH({1987:1988}) = SIA_NH@_FillValue



;*****************
;normalise
;*******************
    
;  SIA_SH = (SIA_SH - avg(SIA_SH))/stddev(SIA_SH)
;  SIA_NH = (SIA_NH - avg(SIA_NH))/stddev(SIA_NH)

;*********************
;plots
;*********************


  wks = gsn_open_wks(oType, diro+filo)


  refyear = 2005

  plot = new(2, graphic)

;total time plot
  res             = True
  res@gsnDraw       = False
  res@gsnFrame      = False
;  res@vpHeightF     = 0.5
;  res@vpWidthF      = 1.
  res@trXMinF       = yrst
  res@trXMaxF       = yren
  res@xyLineThicknessF = 3.

  res@gsnXRefLine = refyear
  res@gsnXRefLineColor = "red"

  

  res@gsnLeftString = "global total SIE"
  res@tiYAxisString = "global SIE (milliun km~S~2~N~)"
  plot(0) = gsn_csm_xy(wks, year, (SIA_SH+SIA_NH), res)

delete([/res@trXMinF, res@trXMaxF/])

  res@xyMarkLineMode     = "Markers"
  res@xyMarker           = 16
  res@trXMaxF            = 13.5; ceil(max((/max(SIA_SH), max(SIA_NH)/)))
  res@trXMinF            = 10.;floor(min((/min(SIA_SH), min(SIA_NH)/)))
  res@trYMinF            = res@trXMinF
  res@trYMaxF            = res@trXMaxF
  res@tiYAxisString      = "Arctic SIE (million km~S~2~N~)"
  res@tiXAxisString      = "Antarctic SIE (million km~S~2~N~)"

  res@gsnXRefLine       := avg(SIA_SH({:refyear-1}))
  res@gsnYRefLine        = avg(SIA_NH({:refyear-1}))
  res@gsnXRefLineColor   = "black"

  res@gsnLeftString = "Arctic vs Antarctic"

  plot(1) = gsn_csm_xy(wks, SIA_SH({:refyear-1}), SIA_NH({:refyear-1}),res)
  res@xyMarkerColor = "red"

  res@gsnXRefLine       := avg(SIA_SH({refyear:}))
  res@gsnYRefLine        = avg(SIA_NH({refyear:}))
  res@gsnXRefLineColor   = res@xyMarkerColor
    res@gsnYRefLineColor =   res@gsnXRefLineColor  

  oplot =  gsn_csm_xy(wks, SIA_SH({refyear:}), SIA_NH({refyear:}),res)
  overlay(plot(1), oplot)


  resP      = True
  resP@gsnMaximize= True
  gsn_panel(wks, plot, (/1,2/), resP)

 
end
