;for given time range, calculate and plot  SIC, atmos anomaly, 
;from monthly obs
;can plot a time-mean or an animation
;ice is plotted as a filled contor with climatological ice edge


begin


  isAni  = False
  isAnom = True   ;plot SIC as anom or abs value
  latmax = -45.
  latmin = -90.

  stdat  = 202305   ;start date (YYYYMM)  
  endat  = 202307 ;;  ;end date (YYYYMM) 
  
  ;climatology start and end dates [2004-2018 is the ARGO climatology)
  clst   = 1979
  clen   = 2008


;Obs data ice data
  filoi = "~/data/NSIDC/raw/SH/monthly/seaice_conc_monthly_sh_197811_202312_v04r00.nc";"+\
;  array_append_record("sh_197811_202312_v04r00.nc",\
;  "icdr_sh_2024"+sprinti("%0.2i",ispan(1,endat%100,1))+"_f18_v02r00.nc", 0)

  oiNam = "cdr_seaice_conc_monthly"


;atmo data
  oaNam = "msl"
  asrc  = "ERA5"
  ;miust be a single file, to avoid short-to-float conversion inconsistencies with ERA data
  dira  = "~/data/REANALYSIS/"+asrc+"/"
  filoa = asrc+"_msl_monthly_SH_"+(/"1979-2023","2024"/)+".nc"




;plot file
  dirp  = "plots/"
  filp  = "NSIDC_sic_"+asrc+"_psl_map_"+stdat+"-"+endat
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



;******************
;read obs atmos
;************************



  fi   := addfiles(dira+filoa,"r") 
  time := cd_calendar(fi[:]->time, -1)  

  dimz   := (/ dimsizes(time),\
  dimsizes(fi[0]->latitude({latmin:latmax})),\
  getfilevardimsizes(fi[0],"longitude")/)
  var     = new(dimz, float)

  st = 0
   


  do i = 0, dimsizes(filoa)-1

    if (dimsizes(getfilevardimsizes(fi[i],oaNam)).eq.4) then
      tmp := dim_max_n_Wrap(short2flt(fi[i]->$oaNam$(:,:,{latmin:latmax},:)),1)
    else
      tmp:= short2flt(fi[i]->$oaNam$(:,{latmin:latmax},:))
    end if


    ntim                  = dimsizes(tmp(:,0,0))
    var(st:st+ntim-1,:,:) = tmp
    st                    = st+ntim

  end do
  delete(tmp)


  var&time := time


  ;get anomaly
  oatmo = var({stdat:endat},:,:)
  nmon  = dimsizes(oatmo(:,0,0))

  mInd =(/ ind(time.eq.clst*100+(stdat%100)), ind(time.eq.clen*100+(stdat%100))/)
 
  do i = 0, nmon-1
    clim         = dim_avg_n(var(mInd(0):mInd(1):12,:,:),0)
    oatmo(i,:,:) = oatmo(i,:,:)  - clim
    mInd         = mInd+1
  end do

  delete([/var,clim/])


  if (oatmo@units.eq."Pa") then
    oatmo  = oatmo * 1.e-2  ;convert to hPa
  end if

  

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

  if (isAnom) then
    gsn_define_colormap(wks,  "cmp_b2r")
    gsn_reverse_colormap(wks)
    res@cnMaxLevelValF         = 1.
    res@cnMinLevelValF         = -res@cnMaxLevelValF
    res@cnLevelSpacingF        = res@cnMaxLevelValF/10

 
  else

    res@cnMaxLevelValF         = .95
    res@cnMinLevelValF         = 0.15
    res@cnLevelSpacingF        = 0.05
    cmap = read_colormap_file("WhiteBlue")    ; 256 x 4 float array
    res@cnFillPalette = cmap(:161:-1,:) 
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
  ires@cnLineThicknessF      = 4.
  ires@cnLevelSelectionMode  = "ExplicitLevels"
  ires@cnLevels              = 0.15
  ires@gsnLeftString         = ""
  ires@gsnRightString        = ""



;atmo plot resurces
  ores                      = True
  ores@cnInfoLabelOn        = False
  res@cnLineLabelsOn        = False
  ores@gsnDraw              = False
  ores@gsnFrame             = False
  ores@gsnLeftString        = ""
  ores@gsnRightString       = ""
;  ores@cnLevelSelectionMode = "ManualLevels"
;  ores@cnMinLevelValF       = -10.
;  ores@cnMaxLevelValF       = -ores@cnMinLevelValF
;  ores@cnLevelSpacingF      = 2.

  ores@gsnContourZeroLineThicknessF = 0.
  ores@gsnContourNegLineDashPattern = 2
  ores@cnLineThicknessF            = 3.5
  ores@cnLineLabelsOn              = False


;;;;;;;;;;;;;;;
  function overlay_plot(sic[3][*][*]:numeric, atmo[*][*]:numeric)
  local plt,i,iplot, aplot

  begin

     ;SIC 
    plt = gsn_csm_contour_map_polar(wks,sic(2,:,:),res)
    
    ;climatological sea ice edge
    ires@cnLineColor = "magenta"
    iplot            = gsn_csm_contour(wks,sic(0,:,:) ,ires)
    overlay(plt,iplot)

    ;month sea ice edge
    ires@cnLineColor = "green4"
    iplot            = gsn_csm_contour(wks,sic(1,:,:) ,ires)
    overlay(plt,iplot)

    ;atmo
    aplot = gsn_csm_contour(wks,atmo , ores)
    overlay(plt,aplot)

    return(plt)
  end  
 ;;;;;;;;;;;;; 



 

  if (isAni) then

    do mm = 0, nmon-1
      res@tiMainString = oatmo&time(mm)

      plot = overlay_plot(osic(:,mm,:,:),oatmo(mm,:,:))  ;obs
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
    ListAppend(oList, dim_avg_n_Wrap(oatmo,0))
    
    res@tiMainString = stdat+" - "+endat
    plot = overlay_plot(oList[0],oList[1])  ;obs

    draw(plot)
    frame(wks)
  end if


end
