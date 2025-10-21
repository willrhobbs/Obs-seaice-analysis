;for given time range, calculate and plot SST from monthly obs
;can plot a time-mean or an animation
;ice is plotted as a white shade with transparency indicating SIC with climatological ice edge


begin


  isAni  = True
  isAnom = True   ;plot SST as anom or abs value
  latmax = -40.


  stdat  = 202302   ;start date (YYYYMM)  
  endat  = 202508 ;;  ;end date (YYYYMM) 
  
  ;climatology start and end dates [2004-2018 is the ARGO climatology)
  clst   = 1982
  clen   = 2011

  dirb   = "/g/data/gv90/wrh581/"

;Obs data ice data
  isrc = "OSISAF"
  diroi = dirb+isrc+"/"; monthly_sic/"
  filoi = systemfunc("ls "+diroi+"ice_conc_sh_ease2-250_*.nc")

;  filoi = diroi+"sic_pss25_"+array_append_record("197811-202412_v05r00.nc",\
; "2025"+sprinti("%0.2i",ispan(1,endat%100,1))+"_F17_icdr_v03r00.nc" ,0)
  
 oiNam = "ice_conc" ;"cdr_seaice_conc_monthly" ; "cdr_seaice_conc_monthly"


;SST data
  oNam  = "sst"
  osrc  = "OISSTv2"
  diro  = dirb+"/OISST/"
  filo  = "OIv2_SST_"+(/"1982-2025"/)+".nc"


;plot file
  dirp  = "plots/"
  filp  = isrc+"_sic_"+osrc+"_sst_map_"+stdat+"-"+endat
  pType = "png"
  if (isAni) then
    filp  = filp+"_ani"
  end if



;********************
;read obs SIC data 
;*************************



  ;get time array and sort in time order

  fi    = addfiles(filoi,"r")
  time = fi[:]->time
  if (typeof(time).eq."int64") then
    ttmp  = time
    time := toint(ttmp)
    copy_VarAtts(ttmp, time)
    delete(ttmp)
  end if
  time  := cd_calendar(time,-1)




  if (.not.isMonotonic(time)) then
    pq = dim_pqsort(time, 2)
    fi = addfiles(filoi(pq),"r")
    delete(pq)
  end if

  ;read data
  tmp   = fi[:]->$oiNam$ ;(time,ygrid,xgrid)

   if (isbyte(tmp)) then
      tmp  := byte2flt(tmp)
   end if 

;correct units
  if (isrc.eq."OSISAF") then
     tmp := int2flt(tmp) 
     tmp = tmp * 1.e-4
  end if

  tmp&time := time

  ;get climatology and anom
  st   = ind(time.eq.stdat)
  en   = ind(time.eq.endat)
  nmon = en-st+1

  dimz          = dimsizes(tmp(0,:,:))
  osic          = new((/2,nmon,dimz(0),dimz(1)/), typeof(tmp))
  osic(1,:,:,:) = tmp({stdat:endat},:,:)   ;anomaly period
  
;climatology
  mInd = stdat%100

  do i = 0, nmon-1

    ;OSISAF doesn't have a monthly mean for Dec 1987 (sensibly) so get correct time indeces else 
    ;climatology is biased
    clInd := ind(time%100.eq.mInd.and.time/100.ge.clst.and.time/100.le.clen)
    osic(0,i,:,:) = dim_avg_n(tmp(clInd,:,:),0)

    if (mInd.eq.12) then
      mInd = 1
    else
      mInd         = mInd+1
    end if 
  end do



;add coords
  osic@lat2d  = fi[0]->lat;itude
  osic@lon2d  = fi[0]->lon;gitude


;******************
;read sst
;************************



  fo    = addfiles(diro+filo,"r") 

  tmp   := fo[:]->$oNam$(:,{:latmax},:)
 
  tmp&time := cd_calendar(fo[:]->time, -1)

  sst       = tmp({stdat:endat},:,:)


  ;get anomaly
  if (isAnom) then

    nmon  = dimsizes(sst(:,0,0))
    mInd :=(/ ind(tmp&time.eq.clst*100+(stdat%100)), ind(tmp&time.eq.clen*100+(stdat%100))/)
 
    do i = 0, nmon-1
      clim         = dim_avg_n(tmp(mInd(0):mInd(1):12,:,:),0)
      sst(i,:,:)   = sst(i,:,:)  - clim
      mInd         = mInd+1
    end do

    delete(clim)
  end if

  delete(tmp)



  

;***********************************
;plot set-up
;contour fill for SIC anom
;colored line for climatological ice edge

;************************************



  wks = gsn_open_wks(pType, dirp+filp)




;SST plot resources
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
  res@lbTitleString       = "Monthly SST Anomaly (degree C) from OISST data"   ; <-- THIS is the colorbar label
  res@lbTitleFontHeightF  = 0.015
  res@lbTitlePosition     = "Bottom"
  res@cnInfoLabelOn          = False
  res@cnLineLabelsOn         = "False"
  res@gsnLeftString          = ""

  if (isAnom) then
    res@gsnRightString         = "("+clst+"-"+clen+" climatology)"
    res@cnMaxLevelValF         = 3.
    res@cnMinLevelValF         = -res@cnMaxLevelValF
    res@cnLevelSpacingF        = res@cnMaxLevelValF/10
    res@cnFillPalette          =  "cmocean_balance"
  else

    res@gsnRightString         = ""
    res@cnMaxLevelValF         =  1.
    res@cnMinLevelValF         = -1.5
    res@cnLevelSpacingF        =  0.1
    res@cnFillPalette          =  "cmocean_thermal"
  
  end if


 
;Ice edge plot resources
  ires                       = True
  ires@gsnDraw               = False
  ires@gsnFrame              = False
  ires@trGridType            = "TriangularMesh"
  ires@cnLineLabelsOn        = False
  ;ires@cnLineLabelString = "Climatological Ice Edge" ; DOESNT WORK
  ires@cnInfoLabelOn         = False
  ires@cnLineThicknessF      = 5.
  ires@cnLevelSelectionMode  = "ExplicitLevels"
  ires@cnLevels              = 0.15
  ires@gsnLeftString         = ""
  ires@gsnRightString        = ""

;sic plot resorces 
  cres                       = ires
  cres@cnLinesOn             = False
  cres@cnFillOn              = True
  cres@cnLevels             := (/0.15,0.3,0.4,0.5,0.6,0.7,0.8,0.9/) 
  nlev                       =   dimsizes(cres@cnLevels)
  imap                       = new((/nlev,4/),float)         
  imap(:,:2)                 = 1.          ;set RGB to white for every level
  imap(:,3)                  = fspan(0.,1.,nlev)  
  cres@cnFillPalette         = imap





;;;;;;;;;;;;;;;
  function overlay_plot(sic[2][*][*]:numeric, sst[*][*]:numeric)
  local plt,i,iplot, aplot

  begin

     ;SST 
    plt = gsn_csm_contour_map_polar(wks,sst,res)
    
    ;climatological sea ice edge
    ires@cnLineColor           = "black"
    iplot            = gsn_csm_contour(wks,sic(0,:,:) ,ires)
    overlay(plt,iplot)

    ;month ice edge
    ires@cnLineColor           = "violet"
    ;ires@cnLineLabelString = "Monthly Ice Edge"
    iplot2            = gsn_csm_contour(wks,sic(1,:,:) ,ires)
    overlay(plt,iplot2)



    ;sic atmo
    aplot = gsn_csm_contour(wks,sic(1,:,:), cres)
    overlay(plt,aplot)

    return(plt)
  end  
 ;;;;;;;;;;;;; 



 

  if (isAni) then

    do mm = 0, nmon-1
      res@tiMainString = sst&time(mm) ;OLD VERSION

    ; Build pretty title
 ;   tval  = sst&time(mm)
 ;   year  = tval / 100
  ;  month = tval - year*100
  ;  res@tiMainString = sprinti("Sea Ice & SST maps for %s %04i", month_name(month-1), year) ; THIS LINE DOESNT WORK

    ;--- Draw overlay plot (SST + ice edges) ---
      plot = overlay_plot(osic(:,mm,:,:),sst(mm,:,:))  ;obs
      draw(plot)
      frame(wks)

    end do

    ;convert png files to a gif
    CMD = "convert "+dirp+filp+"*.png "+dirp+filp+".gif"
    system(CMD)
    system("rm "+dirp+filp+"*.png")  ; remove monthly individual maps = clean up 

  else

    oList = NewList("fifo")
    ListAppend(oList, dim_avg_n_Wrap(osic,1))
    ListAppend(oList, dim_avg_n_Wrap(sst,0))
    
    res@tiMainString = stdat+" - "+endat
    plot = overlay_plot(oList[0],oList[1])  ;obs

    draw(plot)
    frame(wks)
  end if


end
