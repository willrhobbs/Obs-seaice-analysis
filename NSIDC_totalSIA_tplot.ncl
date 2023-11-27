load "~/scripts/NCL/ncl_funcs.ncl"

;read Antarctic total SIA, 
;plot climatology and anoms on a seperate plot

begin

  yrst = 1979
  yren = 2024
  clst = yrst  ;climatology end
  clen = 2008  ;climatology start

  month = "ALL" ;month to plot (string). 
            ;For all months, 'ALL', 
            ;otherwies 3-letter season name  month (JAN, FEB)       ;

  diri = "~/data/NSIDC/"
  fili = "NSIDC_SH_totalSIA_monthly_"+(/"v04r00.nc","v02r00.nc"/)
  iNam = "SIA_cdr"

  diro = "plots/"
  filo = "NSIDC_"+iNam+"_monthly_tplot"
  if (month.ne."ALL") then
    filo = str_sub_str(filo,"monthly",month)
  end if
  oType = "ps"



;*******************
;read data
;********************

  yr1 = min((/yrst,clst/))
  yr2= max((/yren,clen/))

  fi = addfiles(diri+fili, "r")


  if (any(monref.eq.month)) then
    ;read individual month
    mInd = ind(monref.eq.month)+1
    SIA  = fi[:]->$iNam$({yr1*100+mInd:yr2*100+mInd:12})

    SIA&time = SIA&time/100 ;convert to year
    SIA     = SIA - avg(SIA({clst:clen}))
    sd      = stddev(SIA({clst:clen}))
    SIA    := SIA({yrst:yren})

  else

    SIA = fi[:]->$iNam$({yr1*100+1:yr2*100+12})


    clim = new(12,typeof(SIA))
    sd   = clim
    
    do m = 1, 12
      clim(m-1) = avg(SIA({clst*100+m:clen*100+m:12})) 
      sd(m-1)   = stddev(SIA({clst*100+m:clen*100+m:12})) 
      SIA({yrst*100+m::12}) = SIA({yrst*100+m::12}) - clim(m-1)
    end do

    SIA := SIA({yrst*100+1:yren*100+12})
  end if

  




;*********************
;plot anoms
;*********************


  wks = gsn_open_wks(oType, diro+filo)

  res                       = True
  res@gsnMaximize           = True
  res@gsnDraw               = True 
  res@gsnFrame              = True
  res@vpHeightF             = 0.5
  res@vpWidthF              = 1.

  res@tiYAxisString         = "SIA (million km~S~2~N~)"
  res@gsnYRefLine           = 0.
  res@xyLineThicknessF      = 3.
  res@gsnBelowYRefLineColor = "orange"
  res@gsnAboveYRefLineColor = "cyan"
 
  res@gsnRightString        = "";"(clim = "+clst+"-"+clen+")"


  res@tiMainString          = ""




  if (month.ne."ALL") then
    aplot = gsn_csm_xy(wks,SIA&time,SIA,res)
  else

    ;for all months, need to get nicer tickmarks

    nyr                 = yren-yrst+1
    res@trXMaxF         = nyr*12
    res@trXMinF         = 1
    res@tmXBMode        = "Explicit"
    res@tmXBValues     := ispan(1, res@trXMaxF,60)
    res@tmXBLabels     := ispan(yrst,yren,5)
    res@tmXBMinorValues:= ispan(1, res@trXMaxF,12)

    aplot = gsn_csm_xy(wks,ispan(1,dimsizes(SIA),1),SIA,res)
  end if




;**************************
;plot standardized anoms
;************************

  res@tiYAxisString           = "~F33~s"
  res@gsnDraw                 = False
  res@gsnFrame                = False


  if (month.ne."ALL") then
    SIA   = SIA/sd
    aplot = gsn_csm_xy(wks,SIA&time,SIA,res)
  else

    do m = 1, 12
      SIA({yrst*100+m::12}) = SIA({yrst*100+m::12})/sd(m-1)
    end do

    aplot = gsn_csm_xy(wks,ispan(1,dimsizes(SIA),1),SIA,res)
  end if

;add lines for st. dev
  lres                   = True
  lres@gsLineThicknessF  = 2.
  lres@gsLineDashPattern = 2

  xpts = (/1,res@trXMaxF,res@trXMaxF,1/)
  ypts = (/-1.96,-1.96,1.96,1.96/)
  line = gsn_add_polyline(wks, aplot, xpts, ypts, lres)

  draw(aplot)
  frame(wks)


end
