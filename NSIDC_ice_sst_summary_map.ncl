;for given time range, calculate and plot SST from monthly obs
;can plot a time-mean or an animation
;ice is plotted as a white shade with transparency indicating SIC with climatological ice edge


begin


  isAni  = True
  isAnom = True   ;plot SST as anom or abs value
  latmax = -40.


  stdat  = 202101   ;start date (YYYYMM)  
  endat  = 202411 ;;  ;end date (YYYYMM) 
  
  ;climatology start and end dates [2004-2018 is the ARGO climatology)
  clst   = 1982
  clen   = 2011

  dirb   = "/g/data/jk72/wrh581/"

;Obs data ice data
 ; filoi = "~/data/NSIDC/raw/SH/monthly/seaice_conc_monthly_"+\  
filoi = dirb+"NSIDC/monthly_sic/seaice_conc_monthly_"+\
  array_append_record("sh_197811_202406_v04r00.nc",\
  "icdr_sh_2024"+sprinti("%0.2i",ispan(7,endat%100,1))+"_f18_v02r00.nc", 0)

  oiNam = "cdr_seaice_conc_monthly"


;SST data
  oNam  = "sst"
  osrc  = "OISSTv2"
  diro  = dirb+"/OISST/"
  filo  = "OIv2_SST_"+(/"1982-2024"/)+".nc"


;plot file
  dirp  = dirb+"/Obs-seaice-analysis/plots/"
  filp  = "NSIDC_sic_"+osrc+"_sst_map_"+stdat+"-"+endat
  pType = "png"
  if (isAni) then
    filp  = filp+"_ani"
  end if



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


  dimz          = dimsizes(tmp(0,:,:))
  osic          = new((/2,nmon,dimz(0),dimz(1)/), typeof(tmp))
  osic(1,:,:,:) = tmp(st:en,:,:)   ;anomaly period
  
;climatology
  mInd =(/ ind(time.eq.clst*100+(stdat%100)), ind(time.eq.clen*100+(stdat%100))/)

  do i = 0, nmon-1
    osic(0,i,:,:) = dim_avg_n(tmp(mInd(0):mInd(1):12,:,:),0)
    mInd          = mInd+1 
  end do


;add coords
  osic@lat2d  = fi[0]->latitude
  osic@lon2d  = fi[0]->longitude



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
    mInd =(/ ind(tmp&time.eq.clst*100+(stdat%100)), ind(tmp&time.eq.clen*100+(stdat%100))/)
 
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
      res@tiMainString = sst&time(mm)

      plot = overlay_plot(osic(:,mm,:,:),sst(mm,:,:))  ;obs
      draw(plot)
      frame(wks)

    end do

    ;convert png files to a gif
    CMD = "convert "+dirp+filp+"*.png "+dirp+filp+".gif"
    system(CMD)
    system("rm "+dirp+filp+"*.png")  ;clean up 

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
