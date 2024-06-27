load "~/scripts/NCL/ncl_funcs.ncl"




begin

  yrst = 1979
  yren = 2023


  seas = "MJJ" ;season name 


  sNam = "Atlantic"
  xmin = 270.
  xmax = 90.





  diri = "~/data/NSIDC/"
  fili = "NSIDC_SH_icearea_"+(/"v04r00_197811-202312.nc",\
         "v02r00_202401-202405.nc"/)
  iNam = "SIA"

  diro  = "~/plots/ANTARCTIC/NSIDC/"
  filo  = "NSIDC_"+iNam+"_"+sNam+"_"+seas
  oType = "png"



;*******************
;read data
;********************



  fi = addfiles(diri+fili, "r")
  SIA = fi[:]->$iNam$({yrst*100+1:yren*100+12},:)



  SIA = SIA*1.e-6

 ;correct missing values
  SIA({198712:198801},:) = SIA@_FillValue


   printVarSummary(SIA)
;sector integral

  if (xmin.lt.xmax) then
    SIAsec = dim_sum_Wrap(SIA(:,{xmin:xmax}))
  else
    SIAsec = dim_sum_Wrap(SIA(:,{xmin:}))
    SIAsec = SIAsec + (/ dim_sum(SIA(:,{:xmax}))/)
  end if

;get seasonal mean
  SIA := month_to_seasonN(SIAsec,seas)




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
  res@gsnYRefLine          = avg(SIA)


  res@xyMonoDashPattern    = True
  res@xyLineThicknessF     = 4.

  res@trXMinF              = yrst
  res@trXMaxF              = yren+1

  res@tiMainString = seas+" SIA ("+xmin+"-"+xmax+"~S~o~N~E)"
  plot = gsn_csm_xy(wks, ispan(yrst,yren,1), SIA,  res)



end
