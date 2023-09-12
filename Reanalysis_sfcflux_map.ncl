


;map plot of seasonal-mean sfc fluxes for given year (inc anoms)

begin

  yr     = 2023
  clst   = 1979  ;climatology start date
  clen   = 2008  ;climatology end
  seas   = [/"JJA",6,3/] ;season,first month index (1-12),number of months in season)
 

  latmin = -90.
  latmax = -50.


  src    = "ERA5"
  diri   = "~/data/REANALYSIS/"+src+"/"
  iNam   = (/"slhf","sshf","str","ssr","net"/)
  sign   = (/1,1,1,1,1/)  ;sign convention (to ensure +ve is into ocean

  fili   = src+"_sfcflux_monthly_SH_"+(/"1979-2022","2023"/)+".nc"

  
  dirp   = "plots/"
  filp   = src+"_sfcflux_"+seas[0]+"_"+yr+"_map"
  pType  = "ps"





;***************
;read data
;*****************

  fi           = addfiles(diri+fili,"r")
  time         = cd_calendar(fi[:]->time,-1)
  year         = get_unique_values(time/100)
  

  lat          = fi[0]->latitude({latmin:latmax})
  lon          = fi[0]->longitude

  dimz         = (/5,dimsizes(year),dimsizes(lat), dimsizes(lon)/)
  var          = new(dimz, float)
  var!1        = "year"
  var&year     = year
  var!2        = "lat"
  var&lat      = lat
  var!3        = "lon"
  var&lon      = lon


  do i = 0, dimsizes(fili)-1

    ttmp := cd_calendar(fi[i]->time,-1)


    do v = 0, 3
      if (any(getfilevarnames(fi[i]).eq."expver")) then
        tmp  := dim_max_n_Wrap(fi[i]->$iNam(v)$(:,:,{latmin:latmax},:),1)
      else
        tmp  := fi[i]->$iNam(v)$(:,{latmin:latmax},:)
      end if
    
      tmp     := short2flt(tmp)
 
      do y = ttmp(0)/100,max(ttmp)/100
        st             = ind(ttmp.eq.y*100+seas[1])
        var(v,{y},:,:) = dim_avg_n(tmp(st:st+seas[2]-1,:,:),0)
      end do
    end do
  end do


  delete(tmp)

;sign convention
  do v = 0, 3
    var(v,:,:,:) = var(v,:,:,:) * sign(v)
    iNam(v) = fi[0]->$iNam(v)$@long_name ;get a more meaningful var name

  end do



;net
  var(4,:,:,:) = dim_sum_n(var(:3,:,:,:),0)


;convert from J/day/m^2 to W/m^2
  var  = var/(60.^2 * 24)   





;*******************
;map plots
;**********************

  wks = gsn_open_wks(pType, dirp+filp)

  mres                      = True
  mres@gsnDraw              = False
  mres@gsnFrame             = False
  mres@gsnPolar             = "SH"
  mres@mpMaxLatF            = latmax

  mres@cnLinesOn            = False
  mres@cnFillOn             = True
  mres@cnFillPalette        = "BlueWhiteOrangeRed"
  mres@cnFillMode           = "RasterFill"
  mres@cnLevelSelectionMode = "ManualLevels"
  mres@cnMaxLevelValF       = 120.
  mres@cnMinLevelValF       = -mres@cnMaxLevelValF
  mres@cnLevelSpacingF      = 10.
  mres@lbLabelBarOn         = False
  mres@cnInfoLabelOn        = False

  mres@gsnCenterStringFontHeightF = 0.03
  iNam = str_sub_str(iNam,"Surface ","")

  map = new(5, graphic)
  
  do i = 0, 4
    mres@gsnCenterString = iNam(i)
    map(i) = gsn_csm_contour_map_polar(wks,var(i,{yr},:,:),mres)  
  end do



;panel
  resP                   = True
  resP@gsnMaximize       = True
  resP@txString          = yr
  resP@gsnPanelLabelBar  = True
;  resP@lbOrientation     = "Vertical"
  resP@lbBoxEndCapStyle  = "TriangleBothEnds"
  
  resP@gsnPanelRowSpec   = True

  gsn_panel(wks, map, (/2,3/), resP)




;*********************
;calc climatology and anoms
;**************************

  clim = dim_avg_n_Wrap(var(:,{clst:clen},:,:),1)
  var  = var - conform(var,clim,(/0,2,3/))


;climatology plot
  do i = 0, 4
    mres@gsnCenterString = iNam(i)
    map(i) = gsn_csm_contour_map_polar(wks,clim(i,:,:),mres) 
  end do
  resP@txString            := ""+clst(0)+"-"+clen+" mean"
  gsn_panel(wks, map, (/2,3/), resP)



;****************
;anomaly plot
;**************



  mres@cnMaxLevelValF       = 40.
  mres@cnMinLevelValF       = -mres@cnMaxLevelValF
  mres@cnLevelSpacingF      = 2.
  resP@txString            := yr+" anom"

  do i = 0, 4
    mres@gsnCenterString = iNam(i)
    map(i) = gsn_csm_contour_map_polar(wks,var(i,{yr},:,:),mres) 
  end do

  gsn_panel(wks, map, (/2,3/), resP)


end
