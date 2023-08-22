load "~/scripts/NCL/ncl_funcs.ncl"
load "~/scripts/NCL/plot_funcs.ncl"
;for given time range, calculate and plot climatological mean SIC,
;from monthly obs



begin



  latmax = -53.

  clst   = 1985
  clen   = 2015
  month  = (/"FEB","SEP"/)


;Obs data
  filoi = "~/data/NSIDC/raw/SH/monthly/seaice_conc_monthly_sh_197811_202212_v04r00.nc"


  oiNam = "cdr_seaice_conc_monthly"



;plot file
  dirp  = "~/plots/ANTARCTIC/NSIDC/"
  filp  = "NSIDC_cdr_climatology_map"
  pType = "ps"






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
  st   = ind(time.eq.clst*100+1)
  en   = ind(time.eq.clen*100+12)

  sic   = fi[:]->$oiNam$(st:en,:,:) ;(time,ygrid,xgrid)
  sic  := byte2flt(sic)


;get climatology
  nmon = dimsizes(month)
  clim = sic(:nmon-1,:,:)

  do m = 0, nmon-1
    mInd = ind(monref.eq.month(m))
    clim(m,:,:) = dim_avg_n(sic(mInd::12,:,:),0)
  end do

  ;mask ice-free areas
  msk = where(dim_stddev_n(sic,0).eq.0., 0, 1)
  clim  = mask(clim,msk,1)
 
  delete([/sic,msk/])
    
   
;add coords
  clim@lat2d  = fi[0]->latitude
  clim@lon2d  = fi[0]->longitude
  clim!1      = "time"
 
 

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
  res@lbLabelBarOn           = False
  res@cnInfoLabelOn          = False
  res@cnLineLabelsOn         = False



  res@cnMaxLevelValF         = .95
  res@cnMinLevelValF         = 0.15
  res@cnLevelSpacingF        = 0.05
  cmap = read_colormap_file("WhiteBlue")    ; 256 x 4 float array
  res@cnFillPalette          = cmap(:161:-1,:)
  res@mpOceanFillColor       =  cmap(161,:)  ;set all ocean (even missing) to blue
  res@mpFillDrawOrder        = "PreDraw"     ;draw background foirst

  res@gsnRightString          = ""
  res@gsnLeftStringFontHeightF = 0.03


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
  ires@cnLineColor            = "magenta"

  plot = new(nmon, graphic)

  do i = 0, nmon-1
    
    res@gsnLeftString = pref(i)+month(i)


   
    plot(i) = gsn_csm_contour_map_polar(wks,clim(i,:,:),res)
    
    ;climatological sea ice edge

    iplot            = gsn_csm_contour(wks,clim(i,:,:) ,ires)
    overlay(plot(i),iplot)

  end do


;panel

  resP                = True
  resP@gsnMaximize    = True
  resP@gsnPanelLabelBar = True
  resP@txString         = "("+clst+"-"+clen+")"

  gsn_panel(wks, plot, (/1,2/), resP)


end
