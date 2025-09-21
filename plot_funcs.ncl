load "ncl_funcs.ncl"

;colormap_test                   - plots colormap for given cmap and levels
;conf_bivariate(x, y,sig,Npts)         - generate points of an ellipse for a bivariate conf. interval  

;panelRowCol(nplot)      -         finds efficient rows/columns for panel plotting
;ellipse(x,y,xpos,ypos,Npts)        - return x/y points along an ellipse
;scatterplot(wks, x,y, title, line, lcol) - ;high level function that produces scatterplot of bivariate data
;get_font(wks)           -return function code of default font for a defined workspace



;************************************
;generic alphabet for panel captions
;*************************************

 pref = (/"a) ","b) ","c) ","d) ","e) ", "f) ","g) ","h) ","i) ","j) ","k) ","l) ","m) ","n) ","o) ","p) ","q) ","r) ","s) ","t) ","u) ","v) ","x) ", "y) ","z) "/)



 undef("monref")
 monref = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

 ;season names (for plotting and indexing)
 undef("seasref")
 seasref = (/"DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ"/)


;*********************
;colormap_test 
;*******************

;quickly plots a colorpalette with defined contour levels as labels

  undef("colormap_test")

  procedure colormap_test(cmap[1]:string, levels[*]:numeric)

  local nlev, opt, colmap

  begin

    nlev             = dimsizes(levels)
    opt              = True
    opt@LabelStrings = levels

    colmap           = span_color_rgba(cmap,nlev)


    wks = gsn_open_wks("X11","colours")
    draw_color_palette(wks,colmap, opt)

  end

;**********************
;conf_bivariate
;****************************

;generate points to plot ellipse of bivariate conf. interval
;distribution is assumed to be Gaussian

;x - 1d x data
;y - 1d y dta, same size as x
;sig - required conf. interval (0 < sig <1)
;Npts - number of coordinate pairs to generate (min 20 is recommended)
;output is dimensioned (/2,Npts/), where x coords are (0,:), y coords are(1,:)


  undef("conf_bivariate")

  function conf_bivariate(x[*]:numeric,y[*]:numeric,sig[1]:numeric,Npts[1]:integer)

  local chi,  EOF, x0, y0, xax, yax, phi, pi, theta, out

  begin
  
    ;check sizes
    if (dimsizes(x).ne.dimsizes(y)) then
      print("conf_bivariate: x and y must be the same size")
      exit
    end if

    if (sig.ge.1.or.sig.le.0.) then
      print("conf_bivariate: significance level must be between 0 1nd 1")
      exit
    end if


    ;get chi-squared value fro sig-level
    chi = cdfchi_x(sig, 2)

    ;origin of ellipse
    x0     = avg(x)
    y0     = avg(y)


    ;major axis lengths
    EOF      = eofunc((/x,y/), 2, False)  ;get eigenvalues/vectors
    
    xax      = sqrt(chi*EOF@eval(0))  ;x-axis radius
    yax      = sqrt(chi*EOF@eval(1))  ;y-axis radius


    ;angle of rotation 
    ;(if x, y are non-covariant this reduces to 0)
    phi      = atan2(EOF(0,1), EOF(0,0))
    
    ;calculate points using following formula
    ;x = x0 + [ xax.cos(theta).cos(phi) ] +  [ yax.sin(theta).-sin(phi) ]
    ;y = y0 + [ xax.cos(theta).sin(phi) ] + [ yax.sin(theta).cos(phi) ]
    

    out      = new((/2,Npts/),float)   ;output array
    pi       = 3.14159
    theta    =  fspan(0.,2.*pi,Npts)   ;radians around circle (2pi)
    
  
  
   ;x-coords
    out(0,:) = x0 + (xax * cos(theta) * cos(phi)) + (yax * sin(theta) * -sin(phi))
   
   ;y-coords
    out(1,:) = y0 + (xax * cos(theta) * sin(phi)) + (yax * sin(theta) * cos(phi))

     
    return(out)


  end

;**************************
;panelRowCol
;****************************
 
 undef("panelRowCol")
 function panelRowCol(nplot[1]:integer,orientation[1]:string)
 
;calc. efficient rows/columns for panel plotting
;nplot - scalar integer for number of plots on panel
;orientation - "landscape" or "portrait"

 local nmajor, nminor

 begin

   if (.not.any(orientation.eq.(/"landscape","portrait"/))) then
     print("error: funcion panelRowCol()")
     print("argument 1 must be 'landscape' or 'portrait'")
     exit
   end if

    ;number of columns
    nminor = toint(sqrt(nplot))   
 
    if (nplot%nminor.gt.0) then
      nmajor = nplot/nminor+1
    else
      nmajor = nplot/nminor
    end if

    if (orientation.eq."portrait") then
      return((/nmajor,nminor/))
    else
      return((/nminor, nmajor/))
    end if


 end


;****************************
;ellipse
;************************

 undef("ellipse")
 
 function ellipse(xax[*]:numeric,yax[*]:numeric,xpos[*]:numeric,ypos[*]:numeric,Npts[1]:integer)

;returns x and y points of an ellipse for plotting
;(at present does not allow for rotation)
;xax  = length along x-axis
;yax  = length along y-axis (must be same size as xax)
;xpos = x-coordinate of origin (must be same size as xax)
;ypos = y-coordinate of origin (must be same size as xax)
;Npt  = scalar for number of coords to return (>20 is recommended)


 local Nyax,Nxax,Nxpos,Nypos,t,out, tmp, pi

 begin


   ;check dimenions
   Nyax = dimsizes(yax)
   Nxax = dimsizes(xax)
   Nxpos = dimsizes(xpos)
   Nypos = dimsizes(ypos)

   if (any(Nyax.ne.(/Nxax,Nxpos,Nypos/))) then
     print("function ellipse: xax,yax,xpos and ypos must have same dimensions")
     exit
   end if


   out      = new((/2,Nyax,Npts/),float)   ;output array
   pi       = 3.14
   t        = conform_dims((/Nyax,Npts/), fspan(0.,2.*pi,Npts),1)

   ;x-coords
   tmp        = conform(out(0,:,:),xax,0)
   out(0,:,:) = (tmp * cos(t)) + conform(tmp,xpos,0)
   
   ;y-coords
   tmp        = conform(tmp,yax,0)
   out(1,:,:) = (tmp * sin(t)) + conform(tmp,ypos,0)

   return(out)
 end



;******************************
;scatterplot
;***********************************
;high level function that produces scatterplot of bivariate data
;plot is not drawn or framed, this m,ust be done explicitly

 undef("scatterplot")
 function scatterplot(wks, x[*]:numeric, y[*]:numeric, title[1]:string, line[*]:graphic, lcol[1]:string)

;line is a graphic variable/array to hold regline data until end of script

 local plot, rcf, regf, res, resl, gval, R2


 begin




   ;x and y same size
   if (dimsizes(x).ne.dimsizes(y)) then
     print("fatal; scatterplot function - x and y must be same dimension")
     exit
   end if



  ;scatterplot resources
  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False


  res@xyMarkLineModes     = "Markers"             ; choose which have markers
  res@xyMarkers           =  16                   ; choose type of marker  

  res@tiYAxisString       = y@long_name
  res@tiXAxisString       = x@long_name
  res@tiMainString        = title

  res@gsnXRefLine         = 0.
  res@gsnYRefLine         = 0.

  res@trXMinF             = floor(min(x))
  res@trXMaxF             = ceil(max(x))

  res@trYMinF             = floor(min(y))
  res@trYMaxF             = ceil(max(y))


  resl                    = True
  resl@gsLineColor        = lcol


   ;regression line data
   rcf                = regline(x, y)
   regf               = (rcf * (/res@trXMinF,res@trXMaxF/) ) + rcf@yintercept
   res@gsnLeftString  =  "slope = "+sprintf("%5.2f",rcf)    ;left string

   ;correlation coeff.
   R2  = ( escorc(x, y) ) ^2
   res@gsnRightString  = "R~S~2~N~ = "+sprintf("%4.2f", R2)    ;left string
   

   gval = num(.not.ismissing(x).and..not.ismissing(y))
   res@gsnCenterString = "N = "+gval

   ;scatterplot
   plot = gsn_csm_xy(wks, x, y, res)

   ;regline
   line = gsn_add_polyline(wks, plot,(/res@trXMinF,res@trXMaxF/) , regf, resl)

   return(plot)

 end 

;*********************
;get_font(wks)
;********************

;return function code of default font for a defined workspace

 
 undef("get_font")
 function get_font(wks) 
 local txid 

 begin 
   
   txid = create "text" textItemClass wks end create 
   getvalues txid 
   "txFont" : font 
   end getvalues 
   delete(txid) 

   return("F"+font) 

  end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
