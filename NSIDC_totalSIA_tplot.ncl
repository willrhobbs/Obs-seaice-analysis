load "ncl_funcs.ncl"

;read Antarctic total SIA, 
;plot climatology and anoms on a seperate plot

begin

  yrst = 1979
  yren = 2025
  clst = 1979  ;climatology end
  clen = 2018  ;climatology start

  month = "ALL" ;month to plot (string). 
            ;For all months, 'ALL', 
            ;'ANN" for annual-mean
            ;otherwies 3-letter season name (ASO) or  month (JAN, FEB)       ;

  diri = "/g/data/gv90/wrh581/NSIDC/"
  fili = "NSIDC_SH_totalSIA_monthly_"+(/"v05r00.nc","v03r00.nc"/);(/"v04r00.nc","v02r00.nc"/)
  iNam = "SIA_cdr"

  diro = "plots/"
  filo = "NSIDC_"+iNam+"_monthly_tplot"
  if (month.ne."ALL") then
    filo = str_sub_str(filo,"monthly",month)
  end if

  oType = "png"



;*******************
;read data
;********************

  yr1 = min((/yrst,clst/))
  yr2= max((/yren,clen/))

  fi = addfiles(diri+fili, "r")
  YYYYMM = fi[:]->YYYYMM

  if (any(monref.eq.month)) then
    ;read individual month
    mInd = ind(monref.eq.month)+1
    st   = ind(YYYYMM.eq.yr1*100+mInd)
    en   = ind(YYYYMM.eq.yr2*100+mInd)

    SIA  = fi[:]->$iNam$(st:en:12)

    SIA&time := YYYYMM(st:en:12)/100 ;convert to year
 ;   SIA       = SIA - avg(SIA({clst:clen}))
    sd        = stddev(SIA({clst:clen}))
    SIA      := SIA({yrst:yren})

  elseif (month.eq."ANN") then

    st       = min(ind(YYYYMM.ge.yr1*100+1))
    en       = max(ind(YYYYMM.le.yr2*100+12))
    tmp      = fi[:]->$iNam$(st:en)
    SIA      = new(yr2-yr1+1,typeof(tmp))
    SIA!0    = "time"
    SIA&time = ispan(yr1,yr2,1)

      do yy = yr1, yr2
         ii := ind(YYYYMM/100.eq.yy)
         if (dimsizes(ii).lt.12) then
            continue
          end if
          SIA({yy}) = avg(tmp(ii))
       end do 
    sd      = stddev(SIA({clst:clen}))

  else
    st       = min(ind(fi[:]->YYYYMM.ge.yr1*100+1))
    en       = max(ind(fi[:]->YYYYMM.le.yr2*100+12))
    SIA      = fi[:]->$iNam$(st:en)
    SIA&time := YYYYMM(st:en) 

      ;convert time array to year
    if (any(seasref.eq.month)) then  ;seasonl mean

    ;get 3-month seasons

        ;pad final year
      enmon  = max(SIA&time)%100
      if (enmon.ne.12) then 
        SIA := array_append_record(SIA, new(12-enmon, typeof(SIA)),0)
        SIA&time := yyyymm_time(yr1, yr2, "integer") 
      end if
      
      SIA     := month_to_season(SIA, month)

      SIA&time = SIA&time/100
  
      ;SIA     = SIA - avg(SIA({clst:clen}))
      sd      = stddev(SIA({clst:clen}))

      SIA      := SIA({yrst:yren})

    else 
      clim = new(12,typeof(SIA))
      sd   = clim
    
      do m = 1, 12
        clim(m-1) = avg(SIA({clst*100+m:clen*100+m:12})) 
        sd(m-1)   = stddev(SIA({clst*100+m:clen*100+m:12})) 
        SIA({yrst*100+m::12}) = SIA({yrst*100+m::12}) - clim(m-1)
      end do
    end if


    SIA := SIA({yrst*100+1:yren*100+12} )
  end if



;print years ranked in order 
;SIA({1988}) = SIA@_FillValue
;pq = dim_pqsort(SIA, 1)
;print(SIA&time(pq))
;exit

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

  res@tiYAxisString         = month+" sea ice area (million km~S~2~N~)";"SIA (million km~S~2~N~)"

  if (month.eq."ALL") then
      res@gsnYRefLine       = 0.
  else
      res@gsnYRefLine       = avg(SIA({clst:clen}))
  end if
  res@gsnYRefLineDashPattern = 1
  res@xyLineThicknessF      = 3.
  res@gsnBelowYRefLineColor = "orange"
  res@gsnAboveYRefLineColor = "cyan"
 
  res@gsnRightString        = "";"(clim = "+clst+"-"+clen+")"
  res@tiMainString          = ""


  if (month.ne."ALL") then
    res@trXMaxF         = yren
    res@trXMinF         = yrst+1
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

;add polylines showing change in variance

/;
  breakyr = 2006
  
  resl                   = True
  resl@gsLineThicknessF  = 3.
  resl@gsLineDashPattern = 12
  line                   = new(4, graphic)

  xpts    = (/res@trXMinF,ind(SIA&time.eq.breakyr*100+12)/)
  ypts    = 1.96 * stddev(SIA({:breakyr*100+12})) 
  line(0) = gsn_add_polyline(wks,aplot,xpts, (/ypts,ypts/), resl)
  line(1) = gsn_add_polyline(wks,aplot,xpts, -(/ypts,ypts/), resl)

  resl@gsLineColor = "red"
  xpts    = (/ind(SIA&time.eq.(breakyr+1)*100+1), res@trXMaxF/)
  ypts    = 1.96 * stddev(SIA({(breakyr+1)*100+1:})) 
  line(2) = gsn_add_polyline(wks,aplot,xpts, (/ypts,ypts/), resl)
  line(3) = gsn_add_polyline(wks,aplot,xpts, -(/ypts,ypts/), resl)
;/ 

  draw(aplot)
  frame(wks)


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

  xpts := (/1,res@trXMaxF,res@trXMaxF,1/)
  ypts := (/-1.96,-1.96,1.96,1.96/)
  line = gsn_add_polyline(wks, aplot, xpts, ypts, lres)

  draw(aplot)
  frame(wks)


end
