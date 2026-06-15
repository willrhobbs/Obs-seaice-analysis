load "ncl_funcs.ncl"




begin

  yrst = 1979
  yren = 2026


  month = "MAY" ;month to plot (string). 
            ;For all months, 'ALL', 
            ;'ANN" for annual-mean
            ;otherwies 3-letter season name (ASO) or  month (JAN, FEB)       ;

;  sNam = (/"Circum","Ross-Am","Weddell-KH","Bell"/)
;  xmin = (/0.,150.,327.,277./)
;  xmax = (/360.,270.,90.,295./)

  sNam = "Bellingshausen";(/"E. hemi"/);"Weddell", /)
  xmin = 277.; (/330./);, 277./)
  xmax = 295. ; (/90./);, 295./)



  diri = "/g/data/gv90/wrh581/NSIDC/"
  fili = "NSIDC_SH_icearea_"+(/"v06r00_197811-202512.nc",\
         "v04r00_202601-202605.nc"/)
  iNam = "SIA"

  diro  = "~/gd/Obs-seaice-analysis/plots/"
  filo  = "NSIDC_"+iNam+"_Ehemi_"+month
  oType = "png"



;*******************
;read data
;********************



  fi = addfiles(diri+fili, "r")
  SIA = fi[:]->$iNam$({yrst*100+1:yren*100+12},:)



  SIA = SIA*1.e-6

 ;correct missing values
 if (yrst.lt.1987.and.yren.ge.1988) then 
   SIA({198712:198801},:) = SIA@_FillValue
 end if

;sector integral
 nsec = dimsizes(sNam)
 SIAsec = SIA(:,:nsec-1) 

 do i =0, nsec-1
   if (xmin(i).lt.xmax(i)) then
     SIAsec(:,i) = dim_sum_Wrap(SIA(:,{xmin(i):xmax(i)}))
   else
     SIAsec(:,i) = dim_sum_Wrap(SIA(:,{xmin(i):}))
     SIAsec(:,i) = SIAsec(:,i) + (/ dim_sum(SIA(:,{:xmax(i)}))/)
   end if
 end do



;seasonal/monthly mean

 if (any(monref.eq.month)) then
    ;read individual month
    mInd = ind(monref.eq.month)+1
    st   = yrst*100+mInd
    en   = yren*100+mInd

    SIAsec   := SIAsec({st:en:12},:)
    SIAsec&time = SIAsec&time/100 ;convert to year  

  elseif (month.eq."ANN") then


    tmp      = SIAsec({yr1*100+1:y2*100+12},:)
    SIAsec      = new(yrst-yren+1,typeof(tmp))
    SIAsec!0    = "time"
    SIAsec&time = ispan(yr1,yr2,1)

      do yy = yrst, yren
         ii := ind(YYYYMM/100.eq.yy)
         if (dimsizes(ii).lt.12) then
            continue
          end if
          SIAsec({yy},:) = dim_avg_n(tmp(ii,:),0)
       end do 
    

  else
    st       = yrst*100+1
    en       = yren*100+12

      ;convert time array to year
    if (any(seasref.eq.month)) then  ;seasonl mean

    ;get 3-month seasons

        ;pad final year
      enmon  = max(SIAsec&time)%100
      if (enmon.ne.12) then 
        SIAsec := array_append_record(SIAsec, new(12-enmon, typeof(SIAsec)),0)
        SIAsec&time := yyyymm_time(yr1, yr2, "integer") 
      end if
      
      SIAsec     := month_to_season(SIAsec, month)

      SIAsec&time = SIAsec&time/100

      SIAsec      := SIAsec({yrst:yren})

    else 
      clim = new((/12,nsec/),typeof(SIAsec))
    
      do m = 1, 12
        clim(m-1,:) = dim_avg_n(SIAsec({clst*100+m:clen*100+m:12}), 0) 
        SIAsec({yrst*100+m::12},:) = SIAsec({yrst*100+m::12},:) - clim(m-1)
      end do
    end if


  end if


SIAsec = 1.e3 *SIAsec
print(SIAsec({2023:2026},0))
print(avg(SIAsec({1979:2018},0)))
print(stddev(SIAsec({1979:2018},0)))


exit
;*******************
;plot
;*************************

  wks = gsn_open_wks(oType, diro+filo)
 
  res                      = True
  res@gsnDraw              = True
  res@gsnFrame             = True
  res@gsnMaximize          = True
  res@vpHeightF            = 0.5
  res@vpWidthF             = 1.


  res@tiYAxisString        = "SIA (10~S~6~N~ km~S~2~N~)"
  if (nsec.eq.1) then
    res@gsnYRefLine          = avg(SIAsec)
  else
    res@gsnYRefLine          = 0.
    SIA = dim_rmvmean(SIA)
  end if


  res@xyMonoDashPattern    = True
  res@xyLineThicknessF     = 4.
  res@xyLineColors         = (/"black", "magenta","red","orange","green"/)
  
  res@pmLegendDisplayMode  = "Always"
  res@xyExplicitLegendLabels       = sNam+" ("+xmin+"-"+xmax+"~S~o~N~E)"

  res@trXMinF              = yrst
  res@trXMaxF              = yren+1

  res@tiMainString = month+" SIA "
  plot = gsn_csm_xy(wks, SIAsec&time, SIAsec(lon|:, time|:),  res)



end
