load "~/scripts/NCL/ncl_funcs.ncl"




begin

  yrst = 1980
  yren = 2024


  seas = "MJJ" ;season name 


;  sNam = (/"Circum","Ross-Am","Weddell-KH","Bell"/)
;  xmin = (/0.,150.,327.,277./)
;  xmax = (/360.,270.,90.,295./)

  sNam = (/"Weddell", "Bellingshausen"/)
  xmin = (/327., 277./)
  xmax = (/10., 295./)



  diri = "/g/data/gv90/wrh581/NSIDC/"
  fili = "NSIDC_SH_icearea_"+(/"v05r00_197811-202412.nc",\
         "v03r00_202501-202505.nc"/)
  iNam = "SIA"

  diro  = "plots/"
  filo  = "NSIDC_"+iNam+"_sector_"+seas
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

;get seasonal mean
 SIA := SIAsec(lon|:,time|::12)

 do i = 0, nsec-1
   SIA(i,:) = (/month_to_seasonN(SIAsec(:,i),seas)/)
 end do




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
    res@gsnYRefLine          = avg(SIA)
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

  res@tiMainString = seas+" SIA "
  plot = gsn_csm_xy(wks, ispan(yrst,yren,1), SIA,  res)



end
