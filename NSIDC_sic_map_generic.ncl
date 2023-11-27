;for given time range, calculate and plot mean SIC, with climatological ice edge 
;from monthly obs and model
;can plot a time-mean or an animation


begin


  isAni  = False
  isAnom = False   ;plot SIC as anom or abs value
  latmax = -54.

  stdat  = 202308   ;start date (YYYYMM)  
  endat  = 202308 ;;  ;end date (YYYYMM) 
  
  ;climatology start and end dates [2004-2018 is the ARGO climatology)
  clst   = 1979
  clen   = 2008


;Obs data
  filoi = "~/data/NSIDC/raw/SH/monthly/seaice_conc_monthly_"+\
  array_append_record("sh_197811_202306_v04r00.nc",\
  "icdr_sh_20230"+ispan(7,9,1)+"_f18_v02r00.nc",0)

  oiNam = "cdr_seaice_conc_monthly"



;plot file
  dirp  = "plots/"
  filp  = "NSIDC_cdr_map_"+stdat+"-"+endat
  pType = "png"
  if (isAni) then
    dirp  = dirp+"ani/"
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
  osic          = new((/3,nmon,dimz(0),dimz(1)/), typeof(tmp))
  osic(1,:,:,:) = tmp(st:en,:,:)   ;anomaly period


;climatology
  mInd =(/ ind(time.eq.clst*100+(stdat%100)), ind(time.eq.clen*100+(stdat%100))/)

  do i = 0, nmon-1
    osic(0,i,:,:) = dim_avg_n(tmp(mInd(0):mInd(1):12,:,:),0)
    mInd          = mInd+1 
  end do

  if (isAnom) then
    osic(2,:,:,:) = osic(1,:,:,:) - osic(0,:,:,:);anomaly
    osic = where(osic.eq.0, osic@_FillValue, osic)
  end if
  
;add coords
  osic@lat2d  = fi[0]->latitude
  osic@lon2d  = fi[0]->longitude
  osic!1      = "time"
  osic&time   = time(st:en)
 
;mask land flags
  osic = where(osic.gt.1., osic@_FillValue, osic)



;***********************************
;plot set-up
;colored lines ofr climatological and actual ice edge

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

  if (isAnom) then
    gsn_define_colormap(wks,  "cmp_b2r")
    gsn_reverse_colormap(wks)
    res@cnMaxLevelValF         = 1.
    res@cnMinLevelValF         = -res@cnMaxLevelValF
    res@cnLevelSpacingF        = res@cnMaxLevelValF/10
    osic(2,:,:,:) = osic(1,:,:,:) - osic(0,:,:,:);anomaly
    osic = where(osic.eq.0, osic@_FillValue, osic)
 
  else

    res@cnMaxLevelValF         = .95
    res@cnMinLevelValF         = 0.15
    res@cnLevelSpacingF        = 0.05
    cmap = read_colormap_file("WhiteBlue")    ; 256 x 4 float array
    res@cnFillPalette = cmap(:161:-1,:) ;cmap(:131:-1,:) 
    osic(2,:,:,:) = osic(1,:,:,:)
    
  end if


  res@gsnRightString         = "("+clst+"-"+clen+" climatology)"
  res@gsnLeftString          = ""
 

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




;;;;;;;;;;;;;;;
  function overlay_plot(sic[3][*][*]:numeric)
  local plt,i,iplot, aplot

  begin

     ;SIC 
    plt = gsn_csm_contour_map_polar(wks,sic(2,:,:),res)
    
    ;climatological sea ice edge
    ires@cnLineColor = "violetred3"
    iplot            = gsn_csm_contour(wks,sic(0,:,:) ,ires)
    overlay(plt,iplot)

    ;month sea ice edge
    ires@cnLineColor = "black"
    iplot            = gsn_csm_contour(wks,sic(1,:,:) ,ires)
   ; overlay(plt,iplot)


    return(plt)
  end  
 ;;;;;;;;;;;;; 


  if (isAni) then

    do mm = 0, nmon-1
      res@tiMainString = osic&time(mm)

      plot = overlay_plot(osic(:,mm,:,:))  ;obs
      draw(plot)
      frame(wks)

    end do

    ;convert png files to a gif
    CMD = "convert "+dirp+filp+"*.png "+dirp+filp+".gif"
    system(CMD)
 ;   system("rm "+dirp+filp+"*.png")  ;clean up 

  else

    oList = NewList("fifo")
    ListAppend(oList, dim_avg_n_Wrap(osic,1))
    
    res@tiMainString = stdat+" - "+endat
    plot = overlay_plot(oList[0])  ;obs

    draw(plot)
    frame(wks)
  end if



end
