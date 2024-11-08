;simplified version of NSIDC_ice_atmo_summary_map.ncl
;for given time range, calculate and plot  SIC, atmos anomaly, 
;from monthly obs;
;good for "ourtreach" talks
;ice is plotted as a filled contor with climatological ice edge


begin

  latmax = -50.
  latmin = -90.

  stdat  = 202306   ;start date (YYYYMM)  
  endat  = 202308 ;;  ;end date (YYYYMM) 
  
  ;climatology start and end dates [2004-2018 is the ARGO climatology)
  clst   = 1979
  clen   = 2008


;Obs data ice data
  filoi = "~/data/NSIDC/raw/SH/monthly/seaice_conc_monthly_"+array_append_record("sh_197811_202406_v04r00.nc",\
  "icdr_sh_2024"+sprinti("%0.2i",ispan(7,endat%100,1))+"_f18_v02r00.nc", 0)
  
  oiNam = "cdr_seaice_conc_monthly"


;atmo data
  oaNam = (/"u10","v10"/)
  asrc  = "ERA5"
  ;miust be a single file, to avoid short-to-float conversion inconsistencies with ERA data
  dira  = "~/data/REANALYSIS/"+asrc+"/"
  filoa = asrc+"_sfcwnd_monthly_SH_"+(/"1979-2023"/)+".nc";,"2024"/)+".nc"




;plot file
  dirp  = "plots/"
  filp  = "NSIDC_sic_"+asrc+"_wnd_map_"+stdat+"-"+endat
  pType = "png"





;********************
;read obs SIC data 
;*************************

  ;get time array and sort in time order
 

  fi    = addfiles(filoi,"r")
  time  = cd_calendar(fi[:]->time,-1)



  if (.not.isMonotonic(time)) then
    pq = dim_pqsort(time, 2)
    fi = addfiles(filoi(pq),"r")
    delete(pq)
  end if

  ;read data
  tmp   = fi[:]->$oiNam$ ;(time,ygrid,xgrid)
  tmp  := byte2flt(tmp)


  ;get climatology and anom
  st   = ind(time.eq.stdat)
  en   = ind(time.eq.endat)
  nmon = en-st+1


  dimz  = dimsizes(tmp(0,:,:))
  osic  = tmp(st:en,:,:)   ;anomaly period
  
;climatology
  mInd =(/ ind(time.eq.clst*100+(stdat%100)), ind(time.eq.clen*100+(stdat%100))/)


  do i = 0, nmon-1
    osic(i,:,:) =  osic(i,:,:) - dim_avg_n(tmp(mInd(0):mInd(1):12,:,:),0)
    mInd          = mInd+1 
  end do

  osic := dim_avg_n_Wrap(osic,0)


;add coords
  osic@lat2d  = fi[0]->latitude
  osic@lon2d  = fi[0]->longitude

  osic = where(osic.eq.0, osic@_FillValue, osic)
  
;******************
;read obs atmos
;************************



  fi   := addfile(dira+filoa,"r") 



  wnd =  NewList("fifo")

  function read_wind(uNam[1]:string)
  local var, oatmo, nmon, mInd, i, clim
  
  begin
      var       = short2flt(fi->$uNam$(:,{latmin:latmax},:))
      ;var&date = var&date/100 
      var&time = cd_calendar(var&time, -1)


  ;get anomaly 
      oatmo = var({stdat:endat},:,:)
      nmon  = dimsizes(oatmo(:,0,0))

      mInd =(/ clst*100+(stdat%100), clen*100+(stdat%100)/)
  
      do i = 0, nmon-1
        clim         = dim_avg_n(var({mInd(0):mInd(1):12},:,:),0)
        oatmo(i,:,:) = oatmo(i,:,:)  - clim
        mInd         = mInd+1
      end do


      return(dim_avg_n_Wrap(oatmo, 0))
   end

   wnd = [/read_wind(oaNam(0)), read_wind(oaNam(1))/]



  

;***********************************
;plot set-up
;contour fill for SIC anom
;colored lines ofr climatological ice edge
;black lines for atmo
;************************************



  wks = gsn_open_wks(pType, dirp+filp)




;SIC plot resources
  res                        = True
  res@gsnDraw                = False
  res@gsnFrame               = False
  res@gsnMaximize            = True
  res@gsnPolar               = "SH"  
  res@mpMaxLatF              = latmax
 
  res@cnLinesOn              = False
  res@cnFillOn               = True
  res@cnFillDrawOrder        = "PreDraw"
  res@cnFillMode             = "RasterFill" ;Raster fill necessary for high-res
  res@cnLevelSelectionMode   = "ManualLevels"
  res@lbLabelBarOn           = True
  res@cnInfoLabelOn          = False
  res@cnLineLabelsOn         = "False"


  gsn_define_colormap(wks,  "cmp_b2r")
  gsn_reverse_colormap(wks)
  res@cnMaxLevelValF         = 1.
  res@cnMinLevelValF         = -res@cnMaxLevelValF
  res@cnLevelSpacingF        = res@cnMaxLevelValF/10
  res@gsnRightString         = "("+clst+"-"+clen+" climatology)"
  res@gsnLeftString          = ""
 


;atmo plot resurces
  ores                      = True
  ores@gsnDraw              = False
  ores@gsnFrame             = False
  ores@gsnMaximize          = True
  ores@gsnLeftString        = ""
  ores@gsnRightString       = ""
  ores@vcGlyphStyle         = "CurlyVector"
  ores@vcLineArrowThicknessF = 2.
  ores@vcMinDistanceF       = 0.03
  ores@vcRefLengthF         = .2
  ores@vcRefMagnitudeF      = 10.  
  ores@vcLineArrowHeadMaxSizeF = 0.018
  ores@vcMinMagnitudeF         = 8.




  ;SIC 
  plt = gsn_csm_contour_map_polar(wks,osic,res)
    


    ;atmo
  aplot = gsn_csm_vector(wks,wnd[0],wnd[1] , ores)
  overlay(plt,aplot)



  draw(plt)
  frame(wks)


end
