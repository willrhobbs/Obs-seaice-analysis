load "ncl_funcs.ncl"
load "plot_funcs.ncl"

;scatterplot of annula-mean Arctic vs Antarctic SIA


begin

  yrst = 1979
  yren = 2024
  clst = 1979  ;climatology end
  clen = 2018  ;climatology start

     ;

  diri = "/g/data/gv90/wrh581/NSIDC/"
  fils = "NSIDC_SH_totalSIA_monthly_"+(/"v05r00.nc","v03r00.nc"/)
  sNam = "SIE_cdr"

  filn = str_sub_str(fils, "SH", "NH")

  diro = "plots/"
  filo = "NSIDC_Arctic_vs_Antarctic"
  oType = "png"



;*******************
;read SH data
;********************


undef("read_data") 

  function read_data(fpath[*]:string, lnam[1]:string)

  local fi, YYYYMM, st, en, tmp, year, SIA

  begin

      fi = addfiles(fpath, "r")

      SIA      = fi[:]->$sNam$
      SIA&time := fi[:]->YYYYMM

      SIA({198712:198801}) = SIA@_FillValue
       SIA@long_name = lnam

      return(SIA)

    end 

 SIA= [/ read_data(diri+fils, "Antarctic"), \
         read_data(diri+filn, "Arctic") /]

  ListAppend(SIA, SIA[0] + SIA[1])  ;total
  SIA[2]@long_name = "Global"
SIA[2]!0 = "time"
SIA[2]&time = SIA[0]&time



;*****************
;annual mean
;********************



  year     = ispan(yrst, yren, 1)
  SIAann      = new((/3,dimsizes(year)/), typeof(SIA[0]))
  SIAann!1    = "year"
  SIAann&year = year
     
  
  do i = 0, 2
      do y = yrst, yren
        SIAann(i,{y}) = avg(SIA[i]({y*100+1:y*100+12}))
      end do
  end do

;******************
;anom
;*******************



  do i = 0, 2
      do mInd = 0, 11
         m = SIA[i]&time(mInd)%100

         clim = avg(SIA[i]({clst*100+m:clen*100+m:12}))
         SIA[i](mInd::12) = SIA[i](mInd::12) - clim
      end do
   end do



;*********************
;Annual mean plots
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
  res@tiYAxisString = "global SIE (million km~S~2~N~)"
  plot(0) = gsn_csm_xy(wks, SIAann&year, SIAann(2,:), res)

delete([/res@trXMinF, res@trXMaxF/])

  res@xyMarkLineMode     = "Markers"
  res@xyMarker           = 16
  res@trXMaxF            = 13.5; ceil(max((/max(SIA_SH), max(SIA_NH)/)))
  res@trXMinF            = 10.;floor(min((/min(SIA_SH), min(SIA_NH)/)))
  res@trYMinF            = res@trXMinF
  res@trYMaxF            = res@trXMaxF
  res@tiYAxisString      = "Arctic SIE (million km~S~2~N~)"
  res@tiXAxisString      = "Antarctic SIE (million km~S~2~N~)"

  res@gsnXRefLine       := avg(SIAann(0,{:refyear-1}))
  res@gsnYRefLine        = avg(SIAann(1,{:refyear-1}))
  res@gsnXRefLineColor   = "black"

  res@gsnLeftString = "Arctic vs Antarctic"

  plot(1) = gsn_csm_xy(wks, SIAann(0,{:refyear-1}), SIAann(1,{:refyear-1}),res)
  res@xyMarkerColor = "red"

  res@gsnXRefLine       := avg(SIAann(0,{refyear:}))
  res@gsnYRefLine        = avg(SIAann(1,{refyear:}))
  res@gsnXRefLineColor   = res@xyMarkerColor
    res@gsnYRefLineColor =   res@gsnXRefLineColor  

  oplot =  gsn_csm_xy(wks, SIAann(0,{refyear:}), SIAann(1,{refyear:}),res)
  overlay(plot(1), oplot)


  resP      = True
  resP@gsnMaximize= True
  resP@txString = "Annual"
  gsn_panel(wks, plot, (/1,2/), resP)


  delete(res)


;*********************
;monthly anom plots
;*********************


  res             = True
  res@gsnDraw       = False
  res@gsnFrame      = False

  res@gsnLeftString = "global total SIE"
  res@tiYAxisString = "global SIE anom (million km~S~2~N~)"
  res@gsnYRefLine   = 0.  
  res@xyLineThicknessF = 3.


;total time plot
  res@tmXBMode      = "Explicit"
  res@trXMinF       = 0
  res@trXMaxF       = dimsizes(SIA[2])-1
  tdum              = ispan(res@trXMinF, res@trXMaxF, 1)
  st = min(ind(SIA[2]&time%100.eq.1))
  res@tmXBValues    = ispan(st, res@trXMaxF, 60)
  res@tmXBLabels    = ispan( SIA[2]&time(st)/100, yren, 5)  
  res@tmXBMinorValues     = ispan(st, res@trXMaxF, 12)
  
   plot(0) = gsn_csm_xy(wks, tdum, SIA[2], res)


  resm               = True  
  resm@gsnDraw       = False
  resm@gsnFrame      = False

  resm@xyMarkLineMode     = "Markers"
  resm@xyMarker           = 16
  resm@trXMaxF            = 3.5 ;ceil(max((/max(SIA[0]), max(SIA[1])/)))
  resm@trXMinF            = -3.5 ;floor(min((/min(SIA[0]), min(SIA[1])/)))
  resm@trYMinF            = resm@trXMinF
  resm@trYMaxF            = resm@trXMaxF
  resm@tiYAxisString      = "Arctic SIE (million km~S~2~N~)"
  resm@tiXAxisString      = "Antarctic SIE (million km~S~2~N~)"

  resm@gsnXRefLine       = 0.
  resm@gsnYRefLine        = 0.
  

  resm@gsnLeftString = "Arctic vs Antarctic"

  plot(1) = gsn_csm_xy(wks, SIA[0], SIA[1],resm)

 ; refyear = 2025
  ii      = ind(SIA[0]&time/100.ge.refyear)
  resm@xyMarkerColor = "red"
   oplot = gsn_csm_xy(wks, SIA[0](ii), SIA[1](ii),resm)
   overlay(plot(1), oplot)

  resP@txString = "Monthly anomalies ("+clst+"-"+clen+" climatology)"
  gsn_panel(wks, plot, (/1,2/), resP)

 
end
