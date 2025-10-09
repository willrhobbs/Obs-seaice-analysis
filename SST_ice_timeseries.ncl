;timeplot of areal-mean monthly sst anomaly for given domain

;--- Load required NCL libraries
load "ncl_funcs.ncl"

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
  diro  = "/g/data/gv90/wrh581/OISST/"
  filo  = "OIv2_SST_1982-2025.nc"

; Sea ice area data
  yrst = 1979
  yren = 2025
  clst_ice = 1979  ;climatology start
  clen_ice = 2008  ;climatology end

;refyear = yren
  refyear = (/2023,2024, 2025/)  ;years to highlight

  diri = "/g/data/gv90/wrh581/NSIDC/"
  fili = "NSIDC_SH_totalSIA_monthly_"+(/"v05r00.nc","v03r00.nc"/)
  iNam = "SIE_cdr"

;output
  dirp  = "plots/"
  filp  = osrc+"_sst_anom_SO_timeseries"
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


;--------------------------------------------
; Compute SST anomalies (retain time dimension)
;--------------------------------------------
climave   = new(12, typeof(sst))
climsd    = new(12, typeof(sst))
sst_anoms = sst   ; make a full copy to preserve dims & coords

do i = 0, 11
  month = sst&time(i)%100
  stdat = clst*100 + month
  endat = clen*100 + month

  climave(i) = avg(sst({stdat:endat:12}))
  climsd(i)  = stddev(sst({stdat:endat:12}))

  sst_anoms(i::12) = sst(i::12) - climave(i)
end do

; ---------- diagnostics: print min/max and summaries ----------
print("===== SST anomaly diagnostics 1 =====")
printVarSummary(sst_anoms)
print("sst_anom min = ")
print(min(sst_anoms))
print("sst_anom max = ")
print(max(sst_anoms))
print("length sst_anom = " + sprinti("%i", dimsizes(sst_anoms)))


;====================================================
; Read and compute SIA anomalies
;====================================================

print("---- Reading SIA ----")

yr1 = min((/yrst,clst_ice/))
yr2 = max((/yren,clen_ice/))
fi  = addfiles(diri+fili, "r")
YYYYMM = fi[:]->YYYYMM
SIA    = fi[:]->$iNam$
SIA&time := YYYYMM
SIA_2keep    = fi[:]->$iNam$
SIA_2keep&time := YYYYMM

; Compute anomaly for each month
do m = 1, 12
  clim = avg(SIA({clst_ice*100+m:clen_ice*100+m:12}))
  SIA({yrst*100+m::12}) = SIA({yrst*100+m::12}) - clim
end do

; Align SIA time range with SST
common_time = ind(SIA&time.ge.min(sst&time) .and. SIA&time.le.max(sst&time))
sie = SIA(common_time)
sie!0 = "time"
sie&time = sst&time

SIA_2plot=SIA_2keep(common_time)
SIA_2plot!0 = "time"
SIA_2plot&time = sst&time


; ---------- diagnostics: print min/max and summaries ----------
print("===== SIA / sie diagnostics =====")
printVarSummary(sie)
print("sie min = ")
print(min(sie))
print("sie max = ")
print(max(sie))
print("length sie = " + sprinti("%i", dimsizes(sie)))
print("SIA min = ")
print(min(SIA_2plot))
print("SIA max = ")
print(max(SIA_2plot))
; -----------------------------------------------------------

;====================================================
; Plot SST anomaly (original colors) and SIA anomaly
;====================================================

print("---- Plot SST + SIA anomalies ----")

pType@wkWidth  = 2500
pType@wkHeight = 2500
wks = gsn_open_wks(pType, dirp+filp)

;--- Convert yyyymm time to decimal years for proper spacing
time = sst&time/100. + mod(sst&time,100)/12.0   ; e.g., 198109 -> 1981 + 9/12 = 1981.75
sie_time = time   ; align SIA overlay

;--- SST base plot
resSST                       = True
;resSST@trXMinF               = min(time)
;resSST@trXMaxF               = max(time)
resSST@gsnDraw               = False
resSST@gsnFrame              = False
resSST@vpHeightF             = 0.5
resSST@vpWidthF              = 0.75
resSST@tiYAxisString         = "SST anomaly (~S~o~N~C)"
resSST@gsnYRefLine           = 0.
resSST@xyLineThicknessF      = 2.5
resSST@gsnBelowYRefLineColor = "deepskyblue"
resSST@gsnAboveYRefLineColor = "orangered"
resSST@xyLineColor           = "black"
resSST@tiMainString          = "Antarctic sea ice area and sea surface temperature anomalies"

; --- Add pale gray grid ---
resSST@tmXMajorGrid = False
resSST@tmYMajorGrid = True
resSST@tmYMinorGrid = True
resSST@tmXMajorGridLineColor = "gray80"
resSST@tmYMajorGridLineColor = "gray80"
resSST@tmYMinorGridLineColor = "gray80"
resSST@tmXMajorGridThicknessF = 1.0
resSST@tmYMajorGridThicknessF = 1.0
resSST@tmYMinorGridThicknessF = 1.0

;--- X-axis labels (years)
resSST@tmXBMode   = "Explicit"
resSST@tmXBValues = ispan(1982, 2025, 5)  ; one label every 5 year
resSST@tmXBLabels = ispan(1982, 2025, 5)
resSST@tmXBMinorValues = ispan(1982, 2025, 1) + 0.5  ; optional: mid-year minor ticks


;--- Plot SST anomalies
plotSST = gsn_csm_xy(wks, time, sst_anoms, resSST)

;--- Add grey shading 
polyres                  = True
polyres@gsFillColor      = "gray85"     ; pale grey fill
polyres@gsEdgeColor      = "gray80"     ; optional thin outline
polyres@gsEdgeThicknessF = 1.0
polyres@gsFillOpacityF   = 0.35          ; makes slightly transparent (if using NCL ≥6.5)

;for Dec 2016 to Feb 2017
xpoly = (/2016.90, 2017.35, 2017.35, 2016.9/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade1 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jan 2018 to Mar 2018
xpoly = (/2018.1, 2018.55, 2018.55, 2018.1/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade2 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jan 2019 to Mar 2019
xpoly = (/2019.1, 2019.55, 2019.55, 2019.1/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade3 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jan 2020 to Mar 2020
xpoly = (/2019.9, 2020.35, 2020.35, 2019.9/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade4 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jul 2022 to oct 2022
xpoly = (/2022.4, 2022.85, 2022.85, 2022.4/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade5 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jul 2023 to oct 2023
xpoly = (/2023.4, 2023.85, 2023.85, 2023.4/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade6 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)

;for Jul 2024 to oct 2024
xpoly = (/2024.4, 2024.9, 2024.9, 2024.4/)  ; x coordinates for Dec–Feb
ypoly = (/-0.4, -0.4, 0.8, 0.8/)  ; full vertical span
extreme_shade7 = gsn_add_polygon(wks, plotSST, xpoly, ypoly, polyres)


;--- Overlay SIA anomalies
resSIE                  = True
resSIE@gsnDraw          = False
resSIE@gsnFrame         = False
resSIE@xyLineThicknessF = 8
resSIE@xyLineColor      = "black"
resSIE@tiYAxisString    = "SIA anomaly (10^6 km~S~2~N~)"
resSIE@gsnRightAxis     = True

;--- scale SIA to fit on SST axis if needed
plotSIE = gsn_csm_xy(wks, sie_time, sie/10 + 0.5, resSIE)  ; for SIA anomaly
;plotSIE = gsn_csm_xy(wks, sie_time, SIA_2plot/20, resSIE)  ; for SIA 

overlay(plotSST, plotSIE)

;--- Now add red overlay for SIA anomalies from Jan 2016 onwards
; Create logical mask for time >= 2016.0
idx_2016 = ind(sie_time.ge.2016.0)

if (.not.any(ismissing(idx_2016))) then
  resSIE_red                  = True
  resSIE_red@gsnDraw          = False
  resSIE_red@gsnFrame         = False
  resSIE_red@xyLineThicknessF = 2.0
  resSIE_red@xyLineColor      = "red"
  resSIE_red@gsnRightAxis     = True

  plotSIE_red = gsn_csm_xy(wks, sie_time(idx_2016), sie(idx_2016)/10 + 0.5, resSIE_red)
  overlay(plotSST, plotSIE_red)
end if

;--- Draw vertical dashed grey line at January 2016
lnres                     = True
lnres@gsLineColor         = "gray40"    ; darker grey line
lnres@gsLineThicknessF    = 10         ; thickness
lnres@gsLineDashPattern   = 2           ; dashed line style

xline = (/2016.0, 2016.0/)              ; X coordinates
yline = (/-0.4, 0.8/) ; Y range of the line

dum_line = gsn_add_polyline(wks, plotSST, xline, yline, lnres)

draw(plotSST)

;--- Add text labels using normalized coordinates (independent of data space)
txres               = True
txres@txFontHeightF = 0.02
txres@txFontColor   = "black"
txres@txJust        = "CenterCenter"

; Coordinates in normalized device space (0–1) NOT WORKING!
;dum1 = gsn_text_ndc(wks, "Sea Surface Temperature", 0.5, 0.3, txres)
;dum2 = gsn_add_text(wks, plotSST, "Sea Ice Extent", 1997, 0.3, txres)



frame(wks)



;====================================================
; Standardized SST anomalies (kept separate)
;====================================================

;print("---- Plot standardized SST anomalies ----")

;sst_std = sst_anoms    ; copy again to avoid overwrite

;do m = 0, 11
;  sst_std(m::12) = sst_std(m::12) / climsd(m)
;end do

;resStd                   = True
;resStd@gsnDraw           = False
;resStd@gsnFrame          = False
;resStd@tiYAxisString     = "~F33~s"
;resStd@xyLineThicknessF  = 2.0
;resStd@xyLineColor       = "black"
;resStd@gsnYRefLine       = 0.
;resStd@gsnBelowYRefLineColor = "deepskyblue"
;resStd@gsnAboveYRefLineColor = "orangered"

;aplot = gsn_csm_xy(wks, sst_std&time, sst_std, resStd)

;; Add ±1.96 lines
;lres                   = True
;lres@gsLineThicknessF  = 2.
;lres@gsLineDashPattern = 2
;xpts = (/min(sst_std&time), max(sst_std&time)/)
;ypts = (/1.96, 1.96/)
;dumU = gsn_add_polyline(wks, aplot, xpts, ypts, lres)
;ypts = (/-1.96, -1.96/)
;dumL = gsn_add_polyline(wks, aplot, xpts, ypts, lres)

;draw(aplot)
;frame(wks)

end