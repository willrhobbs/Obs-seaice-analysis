load "~/scripts/NCL/ncl_funcs.ncl"

;timeplot of areal-mean monthly sst anomaly for given domain

begin

;spatial domain
  latmin = -65.
  latmax = -45.
  lonmin = 0.
  lonmax = 360.

;climatology
  clst = 1982
  clen = 2011

;SST data
  oNam  = "sst"
  osrc  = "OISSTv2"
  diro  = "~/data/OISST/"
  filo  = "OIv2_SST_1982-2024.nc"

;output
  dirp  = "plots/"
  filp  = osrc+"_midlat_sst_anom_SOcean_tplot"
  pType = "png"


;****************
;read SST data
;*****************

  fo    = addfile(diro+filo,"r") 
  sst   = fo->$oNam$(:,{latmin:latmax},{lonmin:lonmax})
 

;areal average
  wgty  = cos(deg2rad(sst&$sst!1$))
  sst  := wgt_areaave(sst,wgty,1.,0)

  sst!0     = "time"
  sst&time  = cd_calendar(fo->time, -1) ;yyyymm time


;get anomaly
  climave   = new(12, typeof(sst))
  climsd    = climave

  do i = 0, 11

    month = sst&time(i)%100

    stdat    = clst*100+month
    endat    = clen*100+month

    climave(i) = avg(sst({stdat:endat:12}))
    climsd(i)  = stddev(sst({stdat:endat:12}))

    sst(i::12) = sst(i::12) - climave(i)

  end do

;********************
;plot anoms
;**********************


  wks = gsn_open_wks(pType, dirp+filp)

  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False

  res@vpHeightF             = 0.5
  res@vpWidthF              = 1.

  res@tiYAxisString         = "SST anomaly (~S~o~N~C)"
  res@gsnYRefLine           = 0.
  res@xyLineThicknessF      = 2.5
  res@gsnBelowYRefLineColor = "deepskyblue"
  res@gsnAboveYRefLineColor = "orangered"
 
  res@tiMainString          = osrc
  res@gsnRightString        = "(clim = "+clst+"-"+clen+")"
  res@gsnLeftString         = lonmin+"-"+lonmax+"~S~o~N~E : "+latmin+"-"+latmax+"~S~o~N~N"

  yrst = sst&time(0)/100
  yren = max(sst&time)/100
  
  nyr                 = yren-yrst+1
  res@trXMaxF         = nyr*12
  res@trXMinF         = 1
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = ispan(1, res@trXMaxF,60)
  res@tmXBLabels      = ispan(yrst,yren,5)
  res@tmXBMinorValues = ispan(1, res@trXMaxF,12)

  aplot = gsn_csm_xy(wks,ispan(1,dimsizes(sst),1),sst,res)


;map plot showing region
  resm = True
  resm@gsnDraw          = False
  resm@gsnFrame         = False
  resm@gsnPolar         = "SH"
  resm@mpMaxLatF        = max((/latmax+10, -45./))
  map = gsn_csm_map_polar(wks, resm)

  resp  = True
  resp@gsLineThicknessF   = 4.

  if (lonmax-lonmin.lt.360.) then
    xpts = (/lonmin, lonmax, lonmax, lonmin,lonmin/)
    ypts = (/latmin,latmin, latmax,latmax,latmin /)
    pgon  = gsn_add_polyline(wks, map,xpts, ypts, resp)
print("nope")
  else
    pgon = [/gsn_add_polyline(wks, map,(/0.,180.,360./), (/latmin,latmin,latmin/), resp), \
    gsn_add_polyline(wks, map,(/0.,180.,360./), (/latmax,latmax,latmax/), resp)/]
  end if

;Panel
    resP          = True
    resP@gsnMaximize    = True
    gsn_panel(wks, (/aplot , map/), (/1,2/),resP)
   

exit

;**************************
;plot standardized anoms
;************************

  res@tiYAxisString           = "~F33~s"
  res@gsnDraw                 = False
  res@gsnFrame                = False




  do m = 0, 11
    sst(m::12) = sst(m::12)/climsd(m)
  end do

  aplot = gsn_csm_xy(wks,ispan(1,dimsizes(sst),1),sst,res)
 

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
