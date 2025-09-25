load "ncl_funcs.ncl"

;read Antarctic total SIA, 
;plot anoms with a single line for each year

begin

  yrst = 1979
  yren = 2025
  clst = 1979  ;climatology end
  clen = 2008  ;climatology start

;refyear = yren

  refyear = (/2023,2024, 2025/)  ;years to highlight

  diri = "/g/data/gv90/wrh581/NSIDC/"
  fili = "NSIDC_SH_totalSIA_monthly_"+(/"v05r00.nc","v03r00.nc"/)
  iNam = "SIE_cdr"

  diro = "plots/"
  filo = "NSIDC_"+iNam+"_monthly_anoms_byyear"

  oType = "png"



;*******************
;read data
;********************

  yr1 = min((/yrst,clst/))
  yr2= max((/yren,clen/))

  fi = addfiles(diri+fili, "r")
  YYYYMM = fi[:]->YYYYMM

 ; st       = ind(fi[:]->YYYYMM.eq.yr1*100+1)
 ; en       = max(ind(fi[:]->YYYYMM.le.yr2*100+12))
 SIA      = fi[:]->$iNam$
  SIA&time:= YYYYMM
 


;calc anomaly
  do m = 1, 12
    clim = avg(SIA({clst*100+m:clen*100+m:12})) 
    SIA({yrst*100+m::12}) = SIA({yrst*100+m::12}) - clim
  end do
    

;reshape
  year        = ispan(yrst,yren,1)
  SIA_ym      = new((/dimsizes(year),12/), typeof(SIA))
  SIA_ym!0    = "year"
  SIA_ym&year = year
  SIA_ym!1    = "month"
  SIA_ym&month = ispan(1,12,1)

  do yr = yrst, yren
   
    tmp    := SIA({yr*100+1:yr*100+12} )
    tmon   := tmp&time%100
    SIA_ym({yr},{tmon}) = (/tmp/)

  end do



;*********************
;plot anoms
;*********************


  wks = gsn_open_wks(oType, diro+filo)


  res                       = True
  res@gsnMaximize           = True
  res@gsnDraw               = False 
  res@gsnFrame              = False
  res@vpHeightF             = 0.5
  res@vpWidthF              = 1.

  res@tiYAxisString         = "Sea ice extent anomaly (million km~S~2~N~)"
  res@gsnYRefLine           = 0.
  res@xyLineThicknessF      = 2.
  res@xyLineColor           = "gray38"
  res@xyMonoDashPattern     = True
 
;  res@gsnRightString        = "(clim = "+clst+"-"+clen+")"
  res@tiMainString          = ""

  res@tmXBMode              = "Explicit"
  res@trXMinF               = 1
  res@trXMaxF               = 12
  res@tmXBValues            = ispan(2,12,2)
  res@tmXBLabels            = monref(1::2)
  res@tmXBMinorValues       = ispan(1,12,2)


  aplot = gsn_csm_xy(wks,ispan(1,12,1),SIA_ym,res)

  res@xyLineThicknessF      = 6.
  res@xyLineColors         := (/"orange","orange3","orangered","orangered4"/)
  res@pmLegendDisplayMode = "Always"
  res@xyExplicitLegendLabels = refyear
  res@pmLegendOrthogonalPosF = -1.23 
  res@pmLegendParallelPosF = .9
res@pmLegendWidthF = 0.18
res@pmLegendHeightF = 0.15
  oplot = gsn_csm_xy(wks,ispan(1,12,1),SIA_ym({refyear},:),res)
  
  overlay(aplot,oplot)

  draw(aplot)
  frame(wks)






end
