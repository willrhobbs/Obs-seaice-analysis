
load "ncl_funcs.ncl"

;timeplot of areal-mean monthly sst anomaly for given domain

begin

;spatial domain
  latmin = -65.
  latmax = -45.
  lonmin = 0.
  lonmax = 360.

;climatology
  clst = 1982
  clen = 2008

  yrst = 1982
  yren = 2024

;SST data
  oNam  = "sst"
  osrc  = "OISSTv2"
  diro  = "/g/data/gv90/wrh581/OISST/"
  filo  = "OIv2_SST_1982-2025.nc"



;refyear = yren
  refyear = (/2023/)  ;years to highlight

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
sst_anoms = sst({yrst*100+1:yren*100+12})   ; make a full copy to preserve dims & coords

do i = 0, 11

  stdat = clst*100 + i+1
  endat = clen*100 + i+1


  climave(i) = avg(sst({stdat:endat:12}))
  climsd(i)  = stddev(sst({stdat:endat:12}))

  sst_anoms(i::12) = sst_anoms(i::12) - climave(i)
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

yr1 = min((/yrst,clst/))
yr2 = max((/yren,clen/))
fi  = addfiles(diri+fili, "r")
YYYYMM = fi[:]->YYYYMM
SIA    = fi[:]->$iNam$
SIA&time := YYYYMM

; Compute anomaly for each month
do m = 1, 12
  clim = avg(SIA({clst*100+m:clen*100+m:12}))
  SIA({yrst*100+m::12}) = SIA({yrst*100+m::12}) - clim
end do



; Align SIA time range with SST
sie = SIA({yrst*100+1:yren*100+12})


; ---------- diagnostics: print min/max and summaries ----------
print("===== SST anomaly diagnostics =====")
printVarSummary(sst_anoms)
print("sst_anom min = ")
print(min(sst_anoms))
print("sst_anom max = ")
print(max(sst_anoms))
print("length sst_anom = " + sprinti("%i", dimsizes(sst_anoms)))

print("===== SIA / sie diagnostics =====")
printVarSummary(sie)
print("sie min = ")
print(min(sie))
print("sie max = ")
print(max(sie))
print("length sie = " + sprinti("%i", dimsizes(sie)))
; -----------------------------------------------------------

;====================================================
; Plot SST anomaly (original colors) and SIA anomaly
;====================================================

print("---- Plot SST + SIA anomalies ----")

pType@wkWidth  = 2500
pType@wkHeight = 2500
wks = gsn_open_wks(pType, dirp+filp)

;--- SST base plot
resSST                       = True
resSST@gsnDraw               = False
resSST@gsnFrame              = False
resSST@gsnMaximize           = True
resSST@vpHeightF             = 0.5
resSST@vpWidthF              = 1.
resSST@tiYAxisString         = "SST anomaly (~S~o~N~C)"
resSST@gsnYRefLine           = 0.
resSST@xyLineThicknessF      = 2.5
resSST@xyLineColor           = "black"          ; main line black
resSST@gsnBelowYRefLineColor = "deepskyblue"    ; fill below
resSST@gsnAboveYRefLineColor = "orangered"      ; fill above
resSST@tiMainString          = "SST and SIE anomaly timeseries"
resSST@tmXBMode              = "Explicit"
tdum                         = ispan(1,dimsizes(sst_anoms),1)
resSST@trXMinF                  = 1
resSST@tmXBValues            = tdum(::60)
resSST@tmXBLabels            = ispan(yrst, yren, 5)
resSST@tmXBMinorValues       = tdum(::12)
resSST@trYMaxF                = 1.
resSST@trYMinF              = -resSST@trYMaxF


;;--- Overlay SIA line
resSIE                  = True
resSIE@gsnDraw          = False
resSIE@gsnFrame         = False
resSIE@xyLineThicknessF = 5.
resSIE@xyLineColor      = "black"
resSIE@tiYAxisString    = "SIA anomaly (10~S~6~N~ km~S~2~N~)"
;resSIE@gsnRightAxis     = True
resSIE@trYMaxF                = 3.
resSIE@trYMinF              = -resSIE@trYMaxF
resSIE@trYReverse = True  ;invert y-axis so that -ve ice anoms line up with warm anoms

plot = gsn_csm_xy2(wks, tdum, sst_anoms, sie, resSST, resSIE)


draw(plot)
frame(wks)

exit

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
