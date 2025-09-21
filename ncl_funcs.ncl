
;pdfnor_p(mn,sigma,N)   - return parametric Normal distribution for input 1d array
;ZonalMean(x:numeric, lat:numeric, lon:numeric) - get zonal-mean of gridded data, regardless of grid type
;month_to_monyear(x)      - converts any-ranked array of monthly data to (/month,year,.../) format

;bracket_ind(x[*],val[1]) - finds the indices of  1-d monotonic array above and below a given value
;NumberString(str[*]:string) - strips any non-nemric values from a string
;fluxCalc(x,dt)          - Calculate any flux (dX/dt) interpolated to time array of input variable
;iseven(x)               - True if input is an even number, false otherwise
;yresCalc(dlat)          - calc. merid. resolution (kms) for a given lat. res
;xresCalc(dlon,lat)      - calc. zonal resolution (kms) for a given lon res and lat
;identity(dim)           - make a squre identity matrix
     ;factorial(n[1]:integer) - factorial (i.e !) function
;NComb(n[1]:integer)     - number of available conbinations 
;deg2rad
;chkgrid  - check that grids of two arrays match

;equiv_size (x:numeric,  mxlag[1]:integer) - sample size corrected for autocorrlation
;LIstDimz(x[*]:list) - returns the dimension sizes of every element of a list
;ice_sector(ice:numeric,xmin[*]:numeric,xmax[*]:numeric) - calculate SIE sector integrals
;read_csv(fnam[1]:string,head[1]:integer,oType[1]:string) - simple function to read a csv to a given output type

;CeilScale(x[1]:numeric) -similar to built-in ceil() function, but works for numbers > 9 or < 1.
;FloorScale(x[1]:numeric) -similar to built-in floor() function, but works for numbers > 9 or < 1.


;***************************
;useful constants
;**************************

;a list of useful constants for use in calcs

 pi = 3.14159
 Re = 6371220.0  ;mean radius of earth (m)
 g  = 9.80665    ;standard gravity
 
;month names (for plotting and indexing)
 undef("monref")
 monref = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

 ;season names (for plotting and indexing)
 undef("seasref")
 seasref = (/"DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ"/)

;********************
;pdfnor_p
;**********************

;returns the parametric Normal distribution
;useful for plotting on distributions


 undef("pdfnor_p")
 function pdfnor_p(mn[1]:numeric, sigma[1]:numeric, bin_width[1]:numeric)

 local Nbin, xdum, p, norm

 begin
   
   Nbin = 100
   xdum = mn + ispan(-Nbin/2,Nbin/2,1) * bin_width
   
   p = cdfnor_p(xdum,conform_dims(Nbin+1,mn,-1), conform_dims(Nbin+1,sigma,-1))
   

   norm            = 100. * (p(1:) - p(:Nbin-1))
   norm@bin_center = xdum(:Nbin-1)

   return(norm)
 end


;****************************************************
;CeilScale
;****************************************************

;similar to built-in ceil() function, but works for numbers > 9 or < 1.

 undef("CeilScale")
 function CeilScale(x[1]:numeric)
  
  local scl, xceil

  begin

    scl = log10(abs(x))
   
    if (scl.lt.0) then
       scl = floor(scl)
    else
      scl = toint(scl)
    end if

    scl  := 10.^-scl
    xceil = ceil(x*scl)/scl

    return( xceil )

  end



;****************************************************
;FloorScale
;****************************************************

;similar to built-in floor() function, but works for numbers > 9 or < 1.

 undef("FloorScale")
  function FloorScale(x[1]:numeric)
  
  local scl, xfloor

  begin

    scl = log10(abs(x))

    if (scl.lt.0) then
      scl = floor(scl)
    else
      scl = toint(scl)
    end if

    scl   := 10.^scl

    xfloor = scl * floor(x/scl)

    return( xfloor )

  end


;****************************************************
;ZonalMean
;(NB work in progress! - need to work with tricky grids!
;*******************************************************

;get zonal-mean of gridded data, regardless of grid type
;x can be any dimension, but rightmost dims must be y, x
;N.B. It's assumed that var is land/bathymetry-masked
;if binning is required, the lat and lon arrays must each have attributes 'vert', 
;representing the lat and lon vertices respectively, dimensioned (/ny,nx,4/)


  undef("ZonalMean")
  function ZonalMean(var:numeric, lat:numeric, lon:numeric)

  local rank, out, tstx, tsty, dlat, latmin, latmax, nlon, olat, lat_bnds, dimz,i,\
        lat1d,var1d,out1d,ii,yy, lat_vert, lon_vert, wgt

  begin


    ;rectilinear or curvilinear
    rank = dimsizes(dimsizes(lon))
    
    ;simplest case (rectilinear)
    if (rank.eq.1) then
      out = dim_avg_Wrap(var)
      return(out)
    
    else

      ;next simplest case - all lats are in meridional lines
      ;(more likely in SH for tripolar grids)
      tsty = dim_stddev(lat)
      tstx = dim_stddev_n(lon, 0)

      if (all(tstx.eq.0).and.all(tsty.eq.0)) then
        out          = dim_avg_Wrap(var)
        rank         = dimsizes(dimsizes(out))
        out!(rank-1) = "lat"
        out&lat     := lat(:,0)
        return(out)        
      else
        
        if (.not.isatt(lat, "vert").or..not.isatt(lon, "vert")) then
          print("function ZonalMean : for binned calculations lat and lon arguments must have 'vert' attributes")
          exit
        end if

        ;approximate meridional output array
        dlat           = avg(abs( lat(1:,:) - lat(:dimsizes(lat(:,0))-2,:)))
        latmin         = min(lat)
        latmax         = max(lat)
        nlat           = 1+ toint((latmax-latmin)/dlat)
        olat           = fspan(latmin,latmax,nlat)       ;output lat array
        olat@long_name = "latitude"
        olat@units     = lat@units
        dlat           = olat(1)-olat(0)                 ;reset dlat
 

        ;output arrays for binning
        dimz           = dimsizes(var)
        rank           = dimsizes(dimz)
        dimz(rank-2)   = nlat
   
   
        ;reshape arrays to deal with unknown rank
        lat1d          = ndtooned(lat)
        if (rank.eq.2) then
          ;dummy first dimension is required
          var1d        = reshape(var, (/1,product(dimsizes(lat))/))
          out1d        = new((/1,nlat/), typeof(var))          
        else
          var1d        = reshape(var, (/ product(dimz(:rank-3)), dimsizes(lat1d)/))
          out1d        = new( (/ product(dimz(:rank-3)),nlat/), typeof(var))
        end if

        ;output lat bounds (2, nlat)
        lat_bnds       = (/olat-0.5*dlat, olat+.5*dlat/) 

        ;vertices
        lat_vert = reshape(lat@vert, (/dimsizes(lat1d),4/))
        lon_vert = reshape(lon@vert, (/dimsizes(lat1d),4/))
        
        ;loop though dimensions
        do yy = 0, nlat-1

          ii := ind(lat1d.gt.lat_bnds(0,yy).and.lat1d.le.lat_bnds(1,yy))
          if (ismissing(ii(0))) then
            continue
          end if

        ;get areal weights

          ;reset lat vertices so they are within the averaging lat. bound
          yvert      := where(lat_vert(ii,:).gt.lat_bnds(1,yy),lat_bnds(1,yy),lat_vert(ii,:))
          yvert       =  where(yvert.lt.lat_bnds(0,yy),lat_bnds(0,yy),yvert)
          
          wgt        := gc_qarea(yvert,lon_vert(ii,:))

          ;weighted means
          out1d(:,yy) = dim_avg_wgt(var1d(:,ii),wgt,1)

        end do



        ;reshape output array and add coords
        out = (/reshape(out1d,   dimz(:rank-2)) /)
        do i = 0, rank-3
          out!(i)      = var!(i)
          out&$var!(i)$ = var&$var!(i)$ 
        end do
        out!(rank-2)   = "lat"
        out&lat        = olat
        
        return(out)

      end if
    end if

  end


;*************************
;month_to_monyear
;************************

;converts any-ranked array of monthly data to (/month,year,.../) format
;rightmost dimension must be time, and a factor of 12
;attributes and non-time coords are retained 

  undef("month_to_monyear")
  function month_to_monyear(x)

  local dimz, nyr, rank,y,i,span,y3d,x1d,st,j,odimz,y,en

  begin

    dimz = dimsizes(x)
    rank = dimsizes(dimz)

    ;check that time array (rigtmost dimension) is a factor of 12)
    if (dimz(rank-1)%12.ne.0.) then
      print("month_to_monyear(): rightmost dimension must be time and a multiple of 12")
      exit
    end if


    nyr  = dimz(rank-1)/12
   

    if (rank.eq.1) then
      ;simple 1-d case
      y        = new((/12,nyr/),typeof(x))
      do i = 0, 11
        y(i,:) = x(i::12)
      end do

    else
      ;multiple dimensions
      ;use spans and 1d array to avoid if statements with different rank arrays
      span = product(dimz(:rank-2))
      y3d  = new((/span,12,nyr/), typeof(x))
      x1d  = ndtooned(x)

      st = 0
      en = (nyr*12)-1

      do j = 0, span-1       
        do i = 0, 11
          y3d(j,i,:) = x1d(st+i:en:12)
        end do
        
        st = st+(nyr*12)
        en = en+(nyr*12)
      end do
    
    
      ;reshape output array
      odimz = array_append_record(dimz(:rank-2),(/12,nyr/),0)
      y     = new(odimz, typeof(y3d))

      y     = reshape(y3d, odimz)

      ;add coords
      copy_VarCoords_1(x, y)
     
    end if

    ;add metadata
    copy_VarAtts(x, y)
    
    ;time coords
    y!(rank-1) = "month"
    y&month    = ispan(1,12,1)
    y!(rank)   = "year"

    return(y)


  end



;************************
;bracket_ind
;***********************

;finds the indices of  1-d monotonic array above and below a given value
 
 undef("bracket_ind")
 function bracket_ind(x[*]:numeric,val[1]:numeric)

 local tst, st, en

 begin

   tst = isMonotonic(x)

   if (tst.eq.0) then
     print("bracket_in(): x-array must be monotonically increasing or decreasing")
     exit
   end if

   if (tst.eq.1) then
     st = max(ind(x.lt.val))
     en = min(ind(x.gt.val))
   else
     st = min(ind(x.lt.val))
     en = max(ind(x.gt.val))
   end if

   return((/st,en/))
 end


;**************************
;NumberString
;*******************

;strip any no-numeric characters from a string (including minus and decimal points!)

 undef("NumberString")
 function NumberString(str[*]:string)

 local schar, sascii, nstring, out, i, ii

 begin

   schar   = stringtochar(str)  ;convert to character
   sascii  = chartoint(schar)   ;convert to ascii codes

   nstring = dimsizes(str)
   out     = new(nstring, string)

   if (nstring.eq.1) then

     ;make sure there are some numbers
     ;48-57 are characters 0-9 in ascii

     if (any(sascii.ge.48.and.sascii.le.57)) then


       ;get indices of valid numeric charcters
       ;45 is '-;, 46 is '.'
       ii = ind(sascii.ge.48.and.sascii.le.57.or.sascii.eq.45.or.sascii.eq.46)  ;indices of numeric characters

       if (.not.ismissing(ii(0))) then
         out = chartostring(integertochar(sascii(ii)))
       end if
   
     end if

   else
   
     do i = 0, nstring-1

       if (.not.any(sascii(i,:).ge.48.and.sascii(i,:).le.57)) then
         continue
       end if

       ii = ind(sascii(i,:).ge.48.and.sascii(i,:).le.57.or.sascii(i,:).eq.45.or.sascii(i,:).eq.46)  ;indices of numeric characters

       if (.not.ismissing(ii(0))) then
         out(i) = chartostring(integertochar(sascii(i,ii)))
       end if
     
       delete(ii)
     end do
     
   end if

   return(out)

 end



;**************************
;ice_sector
;**************************

;calculate SIE sector integrals
 
 undef("ice_sector")
 function ice_sector(ice:numeric,xmin[*]:numeric,xmax[*]:numeric)

 local nsec, dimz, rank, dNams, out, i

 begin

   ;test dimensions
   nsec = dimsizes(xmin)
   if (dimsizes(xmax).ne.nsec) then
     print("Error : ice_sector() - xmin and xmax must be the same size")
     exit
   end if

   dimz  = dimsizes(ice)
   rank  = dimsizes(dimz)
   dNams = getvardims(ice)

   if (all(dNams(rank-1).ne.(/"lon","longitude"/))) then
     print("Error : ice_sector() - rightmost dimension of first argument must be longitude")
     exit
   end if

   dimz(rank-1) = nsec
   out          = new(dimz,typeof(ice))


   if (rank.eq.1) then
     do i = 0, nsec-1
       if (xmin(i).lt.xmax(i)) then
         out(i) = dim_sum(ice({xmin(i):xmax(i)}))
       else
         out(i) = dim_sum(ice({xmin(i):})) + dim_sum(ice({:xmax(i)}))
       end if
     end do
   end if
   
   if (rank.eq.2) then
     do i = 0, nsec-1
       if (xmin(i).lt.xmax(i)) then
         out(:,i) = dim_sum(ice(:,{xmin(i):xmax(i)}))
       else
         out(:,i) = dim_sum(ice(:,{xmin(i):})) + dim_sum(ice(:,{:xmax(i)}))
       end if
     end do
   end if

   if (rank.eq.3) then
     do i = 0, nsec-1
       if (xmin(i).lt.xmax(i)) then
         out(:,:,i) = dim_sum(ice(:,:,{xmin(i):xmax(i)}))
       else
         out(:,:,i) = dim_sum(ice(:,:,{xmin(i):})) + dim_sum(ice(:,:,{:xmax(i)}))
       end if
     end do
   end if
   

   if (rank.eq.4) then
     do i = 0, nsec-1
       if (xmin(i).lt.xmax(i)) then
         out(:,:,:i) = dim_sum(ice(:,:,:,{xmin(i):xmax(i)}))
       else
         out(:,:,:,i) = dim_sum(ice(:,:,:,{xmin(i):})) + dim_sum(ice(:,:,:,{:xmax(i)}))
       end if
     end do
   end if


   return(out)

 end
 

;**************************
;ListDimz
;****************************

;returns the dimension sizes of every element of a list

 undef("ListDimz")
 function ListDimz(x[*]:list)
 
 local nx, N, i
 
 begin

   nx = ListCount(x)
   N  = new(nx,integer)

   do i = 0, nx-1
     N(i) = sum(dimsizes(x[i]))
   end do

   return(N)

 end 

;********************************
;any_multcond
;******************************

;Returns True for all values of x that equal any values of y.


 undef("any_multcond")
 function any_multcond(x[*],y[*])

 local out, i, ii

 begin

   out = new(dimsizes(x),logical)
   out = False
   
   do i = 0, dimsizes(y)-1

     ii = ind(x.eq.y(i))
     if (.not.ismissing(ii(0))) then
       out(ii) = True
     end if
     
     delete(ii)
   end do

   return(out)

 end


;************************************************
;fluxCalc - 
;************************************************
;calculate flux from input variable
;output is on same time array as input
;first and last entries will be missing values


 undef("fluxCalc")
 function fluxCalc(x[*]:numeric, dt[*]:numeric)

 local ntim, Xint, dX, dXdt

 begin

  
   ; input data, interpolated to input time bounds
   ntim        = dimsizes(x) 
   Xint        = new(ntim,typeof(x))   
   Xint(1:)    = 0.5 * (x(1:) + x(:ntim-2) )
  

   ;difference
   dX          = new(ntim,typeof(x))   ; array for flux difference
   dX(:ntim-2) = Xint(1:) - Xint(:ntim-2)
   
   ;flux
   dXdt        = dX/dt

   return(dXdt)


 end 



;***********************************************
;identify even or odd numbers
;***********************************************

 undef("iseven")
 function iseven(x:numeric)

 local xint, y, out

 begin

   xint = tointeger(x)
   y    = xint%2

   out  = where(y.eq.1, False, True)
   return(out)

 end 


;*******************************
 ;convert degrees to radians
;******************************


 undef("deg2rad")
 function deg2rad(x:numeric)

 local out
 
 begin
   out = (pi/180.) * x
   return(out)
 end

;*****************
;yres(dlat)
;******************

 undef("yresCalc")
 function yresCalc(dlat[*]:numeric)

;return meridional resolution (km) from latitude diffs (degrees)

 local dy

 begin

   dy      = (dlat/180.)* pi * Re * 1.e-3
   return(dy)

 end 


;*****************
;xresCalc(dlon,lat)
;******************

 undef("xresCalc")
 function xresCalc(dlon[*]:numeric,lat[*]:numeric)

;return zonal resolution (km) at a given lat. from longitude diffs (degrees)

 local y, dx, nlat,nlon

 begin

   ;check dimensions
   nlat = dimsizes(lat)
   nlon = dimsizes(dlon)
   
   if (nlat.gt.1.and.nlon.gt.1) then  ;if neither is a scalar
     if (nlat.ne.nlon) then
       print("function xresCalc(): either argument must be a scalar, or both arguments must be the same dimension")
       exit
     end if
   end if
   dx      = (dlon/180.)* pi * Re * cos(deg2rad(lat)) *  1.e-3
   return(dx)
 end 






;*******************************
;identity(dim)
;******************************

 undef("identity")
 function identity(dim[1]:integer)

;create a square identity matrix (dim x dim), type integer

 local I, i

 begin

   I = conform_dims((/dim,dim/),0,-1)
   do i = 0, dim-1
     I(i,i) = 1
   end do
   
   return(I)

 end 




;****************************
;read_csv
;***************************

 undef("read_csv")
 function read_csv(fnam[1]:string,head[1]:integer,oType[1]:string)

;simple function to read a csv to a given output type
;fnam : scalar string representing input filepath
;head : integer value of number of lines to skip at head (if no header set to 0)
;oType : scalar of output type

 local dat, ndat, line, ncol,out,i

 begin

   dat    = asciiread(fnam,-1,"string")
   dsplit = str_split_csv(dat(head:),",",0)
 
 
   if (oType.eq."string") then
     out = dsplit
   else if (oType.eq."float") then
     out = stringtofloat(dsplit)
   else if (oType.eq."integer") then
     out = stringtointeger(dsplit)
   else if (oType.eq."double") then
     out = stringtodouble(dsplit)
   end if
   end if
   end if
   end if

   return(out)
   
 end 


;************************
;factoral (!)
;*************************

 undef("factorial")
 function factorial(n[1]:integer)


 begin 

   if (n.gt.16) then
     print("function ncl_funcs.ncl:factorial - n must be less than 17")
     exit
   end if

   if (n .le. 1) then 
     return(1) 
   else 
     return(n * factorial(n-1)) 
   end if 
 end 


;****************************************
;NComb ; number of available combinations
;*****************************************

 undef("NComb")
 function NComb(n[1]:integer)

 local C, tmp, i

 begin

   C = 0
   
   do i = 1, n

     tmp = factorial(n)/(factorial(i) * factorial(n-i))

     C   = C+tmp
   end do


   return(C)
 end




;*************************
;chk_grid
;*************************

 

 undef("chkgrid")
 function chkgrid(x:numeric, y:numeric, dNam[*]:string)


 ;compare grids of input arrays x, y along dimensions specified by 'dNam'
 ;returns a logical expression, 'True' is they match, 'False' otherwise

;x, y must have same coordinate names


 local rank, rankx, ranky, out, c1x, c1y, tst1, c2x, c2y, tst2, c3x, c3y, tst3 

 
 begin


   ;check array sizes
   rank  = dimsizes(dNam)
   rankx = dimsizes(dimsizes(x))
   ranky = dimsizes(dimsizes(y))
   
   if (rankx.lt.rank.or.ranky.lt.rank) then
     print("function chkgrid: input arrays do not match test dimensions")
     exit
   end if


   ;first coord set
   out = False  ;set low expectations !

   c1x = x&$dNam(0)$
   c1y = y&$dNam(0)$

   if (dimsizes(c1x.eq.dimsizes(c1y))) then  ;check that sizes match
     tst1 = c1x-c1y           ;if they do, check that the locations match

     if (all(abs(tst1).eq.0.)) then      
       out = True
     end if
   end if


   ;2nd coord set
   if (rank.gt.1) then

     c2x = x&$dNam(1)$
     c2y = y&$dNam(1)$

     if (dimsizes(c2x.eq.dimsizes(c2y))) then
       tst2 = c2x-c2y
       if (all(abs(tst2).eq.0.)) then
         out = True
       end if
     end if    
   end if
 

   ;3rd coord set
   if (rank.eq.3) then

     c3x = x&$dNam(2)$
     c3y = y&$dNam(2)$

     if (dimsizes(c3x.eq.dimsizes(c3y))) then
       tst3 = c3x-c3y
       if (all(abs(tst3).eq.0.)) then
         out = True
       end if
     end if    
   end if
       

   return(out)

   end 

;------------------------------------
; =====================================
undef("clmDay1D")
function clmDay1D (x[*]:numeric, yyyyddd:integer)   

; -------

local ntim, ndys, days, clmDay, ndy, indx, year_day, nFill
begin


   ntim  = dimsizes(x)
   ndys  = 366                                     ; allow for leap year

   days  = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   clmDay= new(ndys,typeof(x), getFillValue(x) ) ; daily climatology

;
; Compute averages for each sequebtial day of the year. 
; This uses dimension swapping.
;
   do ndy=0,ndys-2                                 ; ndy=0->364 ==> day of year 1->365       
      indx        = ind( days.eq.(ndy+1) )      
      clmDay(ndy) = avg(x(indx))
      delete(indx)                  
   end do
                                                   ; nominal day 366 
                                                   ; ave(31 Dec + 1 Jan)=leap
   clmDay(ndys-1) = (clmDay(0) + clmDay(ndys-2))*0.5

   nFill = num(ismissing(clmDay))
   if (nFill.eq.0) then
       delete(clmDay@_FillValue)
   end if

   clmDay@long_name   = "Daily Climatology"
   if (isatt(x,"long_name")) then
       clmDay@long_name = clmDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
       clmDay@units     = x@units
   end if
   clmDay@information = "Raw daily averages across all years"
   clmDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   clmDay!0           = "year_day"
   clmDay&year_day    =  year_day
   

   delete(clmDay@year_day)                   ; clean up
 
   return (clmDay)
end

; =====================================
undef("smthClmDay1D")
function smthClmDay1D (clmDay[*]:numeric, nHarm:integer)
;
local nFill, dn, cf, clmDaySmth 
begin
  
  if (isatt(clmDay, "_FillValue")) then
      nFill = num(ismissing(clmDay))
      if (nFill.gt.0) then
          print("smthClmDay_TLL: No missing values allowed: ezfftf does not allow")
          print("smthClmDay_TLL: nFill="+nFill)
          exit  
      end if
  end if



  dn = getvardims(clmDay)                   ; get dimension names
  if (dn(0).ne."year_day") then
      print("smthClmDayTLL: Warning: Usually expect year_day to be the dimension name")
  end if


  cf = ezfftf( clmDay )                        

                                            ; remember NCL is 0-based 
                                            ; cf(:,0:nHarm-1) are retained unaltered
  cf(:,nHarm   ) = 0.5*cf(:,nHarm)  ; mini-taper
  cf(:,nHarm+1:) = 0.0                  ; set all higher coef to 0.0


  clmDaySmth = ezfftb( cf, cf@xbar)                 ; reconstructed series
  clmDaySmth@information = "Smoothed daily climatological averages"
  clmDaySmth@smoothing   = "FFT: "+nHarm+" harmonics were retained."

  return(clmDaySmth)
end

; =====================================
undef("calcDayAnom1D")
function calcDayAnom1D (x[*]:numeric, yyyyddd:integer, clmDay[*]:numeric)   



local ntim, ndys, days, xAnom, nt
begin

   ntim  = dimsizes (x)
   ndys  = 366                               

   days  = yyyyddd - (yyyyddd/1000)*1000      ; strip year info [yyyy]

; quick error check

   if (.not.all(days(0:ndys-2).eq.clmDay&year_day(0:ndys-2))) then
       print("calcDayAnomTLL: day mismatch")
       exit
   end if

; loop on every day and compute difference. 

   xAnom = (/ x  /)                                 ; create xAnom
   do nt=0,ntim-1
      xAnom(nt) = x(nt) - clmDay(days(nt)-1)   ; -1 for 0-based subscript
   end do

   if (isatt(x,"long_name")) then
       xAnom@long_name = "Anomalies: "+x@long_name
   else
       xAnom@long_name = "Anomalies from Daily Climatology"
   end if
   if (isatt(x,"units")) then
       xAnom@units = x@units
   end if

;   copy_VarCoords(x, xAnom)

   return(xAnom)
end

;*******************************
;trajectory_rk4
;******************************

;use 4th order runge-kutte method to integrate a particle's displacement

  undef("trajectory_rk4")

  function trajectory_rk4(u:numeric,v:numeric,lat[*]:numeric, lon[*]:numeric, x0[*]:numeric,y0[*]:numeric,dt[1]:numeric,isCyc[1]:logical)

;  u, v   - velocity fields. Must have same dims and rightmost dims must be (lat,lon)
;  lat, lon - self explanatory
;  x0, y0 - initial lon/lat locations: must be 1d and the same size
;  dt     - scalar; time period (s) of integration
; isCyc - cyclic lon array?  

;return value is an arry of size (/2,dimsizes(x0) /) with lat, lon coords
 
  local dimz, rank, dx, dy, uf, vf,utmp, vtmp, xtmp, ytmp, x1, y1, out 


begin



  ;dim. sizes
  dimz = dimsizes(u)
  if (any(dimz.ne.dimsizes(v))) then
    print("function trajectory_rk4: u, v must be the same size")
    exit
  end if
  if (dimsizes(x0).ne.dimsizes(y0)) then
    print("function trajectory_rk4: x0, y0 must be the same size")
    exit
  end if

  rank = dimsizes(dimz)
  if (dimz(rank-2).ne.dimsizes(lat)) then
    print("function trajectory_rk4: lat array does not correspond to u, v array size")
    exit
  end if
  if (dimz(rank-1).ne.dimsizes(lon)) then
    print("function trajectory_rk4: lon array does not correspond to u, v array size")
    exit
  end if


;1 degree lat/lon in m
  dx = 2*pi*Re*cos(deg2rad(y0))/360.
  dy = 2 * pi * Re/360.

;first iteration
  uf = linint2_points(lon, lat, u, isCyc, x0, y0, 1)
  vf = linint2_points(lon, lat, v, isCyc, x0, y0, 1)

;2nd iteration 
  xtmp = x0 + (dt*uf/dx)/2.  ;locns
  ytmp = y0 + (dt*vf/dy)/2.

  utmp = linint2_points(lon, lat, u, isCyc, xtmp, ytmp, 1)
  vtmp = linint2_points(lon, lat, v, isCyc, xtmp, ytmp, 1)

  uf   = uf + (2*utmp)
  vf   = vf + (2*vtmp)



;3rd iteration'
  xtmp = x0 + (dt*utmp/dx)/2.  ;locns
  ytmp = y0 + (dt*vtmp/dy)/2.

  utmp = linint2_points(lon, lat, u, isCyc, xtmp, ytmp, 1)
  vtmp = linint2_points(lon, lat, v, isCyc, xtmp, ytmp, 1)

  uf   = uf + (2*utmp)
  vf   = vf + (2*vtmp)


;4th iteration
  xtmp = x0 + (dt*utmp/dx)  ;locns
  ytmp = y0 + (dt*vtmp/dy)

  utmp = linint2_points(lon, lat, u, isCyc, xtmp, ytmp, 1)
  vtmp = linint2_points(lon, lat, v, isCyc, xtmp, ytmp, 1)

  uf   = uf + utmp
  vf   = vf + vtmp

;wgted avg

  uf = uf/6.
  vf = vf/6.

  x1 = x0 + (dt*uf/dx)
  y1 = y0 + (dt*vf/dy)

  out = (/x1, y1/)

  return(out)

end



;**********************************
;day_to_seasonN
;*******************************

;convert daily data into 3-month seasonal subsets
undef("day_to_seasonN")
 function day_to_seasonN(x[*][*][*]:numeric, seas[*]:string, time[*]:numeric)

;x         - daily data to transform. Time must be rightmost dim. At the moment only 3-dims are allowed. 
;seas      - 3-month seasons of interest
;time      - coord array of x; must be same as rightmost dim of x and have 'units attribute


 local time, tbase, h, m, s, yrst, yren, dimz, rank, nseas, nyr, out, yrs, monRef, edRef, i, j, mind, emon, eday, st, en, tmp1, tmp2

 begin



   ;time coord info


   if (.not.isatt(time, "units") ) then
     print("fatal: ~/ncl_funcs.ncl - day_to_seasonN: time array must have 'units' attribute")
     exit
   end if

     

   tbase   = floattoint(cd_calendar(time(0), 0) )
   h       = tbase(:,3)
   m       = tbase(:,4)
   s       = tbase(:,5)

   yrst    = tbase(:,0)

   tbase   = floattoint(cd_calendar(time(dimsizes(time)-1), 0) )
   yren    = tbase(:,0)


   ;output array
   dimz    = dimsizes(x)
   rank    = dimsizes(dimz)
   
   if (dimz(rank-1).ne.dimsizes(time)) then
     print("fatal: ~/ncl_funcs.ncl - day_to_seasonN: time array must be same size as rightmost dimension of x")
     exit
   end if

   nseas   = dimsizes(seas)
   nyr     = yren-yrst+1
   yrs     = ispan(yrst, yren, 1)

   out     = new( (/nseas,dimz(0), dimz(1), nyr, 92/), typeof(x) )
   
   out!0   = "season"
   out!1   = x!0
   out!2   = x!1
   out!3   = "year"
   out!4   = "day"
 
   ;season x-refs
   monRef  = (/"DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ"/) 
   edRef   = (/28,31,30,31,30,31,31,30,31,30,31,31/)



   do i = 0, nseas-1


     mind = ind(monRef.eq.seas(i))

     emon = mind+2            ;final month of season

     eday = edRef(mind)       ;final day of final month

     ;array of start times for each year
     st   = ut_inv_calendar(yrs, conform(yrs,mind,-1),conform(yrs,1,-1),conform(yrs,h,-1),conform(yrs,m,-1),conform(yrs,s,-1), time@units, 0)

     ;array of finish times for each year
     en   = ut_inv_calendar(yrs, conform(yrs,emon,-1),conform(yrs,eday,-1),conform(yrs,h,-1),conform(yrs,m,-1),conform(yrs,s,-1), time@units, 0)




      ;DJF
     if (mind.eq.0) then  

       st    = ut_inv_calendar(yrs-1, conform(yrs,12,-1),conform(yrs,1,-1),conform(yrs,h,-1),conform(yrs,m,-1),conform(yrs,s,-1), time@units, 0) 
            ;ignore first 'December' 
       st(0) =  ut_inv_calendar(yrs(0),1,1,h,m,s, time@units, 0)   
     end if




    ;NDJ
     if (mind.eq.11) then

       en   = ut_inv_calendar(yrs+1, conform(yrs,1,-1),conform(yrs,31,-1),conform(yrs,h,-1),conform(yrs,m,-1),conform(yrs,s,-1), time@units, 0) 

         ;ignore last 'January'
       en(nyr-1) = ut_inv_calendar(yrs(nyr-1),12,31,h,m,s, time@units, 0)
     end if


     if (isinteger(time)) then

       tmp1 = st
       tmp2 = en
       delete(st)
       delete(en)
       st   = doubletoint(tmp1)
       en   = doubletoint(tmp2)
       delete(tmp1)
       delete(tmp2)
     end if



     do j = 0, nyr-1


       tmp                                  =  (/ x(:,:,{st(j):en(j)}) /)


       out(i,:,:,j,:dimsizes(tmp(0,0,:))-1) =  (/ tmp /)
       delete(tmp)

     end do

     delete(st)  ;remove, just in case...
     delete(en)
   end do


   return(out)

 end 




;********************************
;get_transect
;**********************************

;returns coord indices for a transect between two points on a grid
 undef("get_transect")
 function get_transect(x1[1]:integer, y1[1]:integer, x2[1]:integer, y2[1]:integer)

;x1, y1 = x (lon) and y (lat) coord indices of second point
;x2, y2 = x (lon) and y (lat) coord indices of second point

;x1 MUST be < x2


 local out, x, N, dx, dy, dydx, y

 begin

   if (x1.gt.x2) then
     print("Error in ~/ncl_funcs.ncl function get_transect. x1 > x2")
     exit
  end if
     

   x         = ispan(0, x2-x1, 1)  ; array of x indices
   N         = dimsizes(x)       ; number of points in line


   ;output array

   if (N.eq.1) then  ;both points on same x line
     
     out     = (/ (/y1, y2/), (/x1,x2/) /)
     
   else

     out       = new((/2, N/), integer)

     out(1,:)   = x+x1    ;all x vals     
     out(0,0)   = y1      ;first point (point 1)
     out(0,N-1) = y2      ;last point (point 2)

     
   ;calc differences  and gradient
     dx        = int2flt(x2 - x1)
     dy        = int2flt(y2 - y1)
     dydx      = dy/dx

   ;nearest indices to straight line
     y         = round( x(1:N-2) * dydx, 3)

     out(0,1:N-2) = y + y1
   end if
     

   return(out)
   


 end 




;**************************************
;rmAnnCycle1D_smth
;*******************************************************************
; 
; remove smoothed annual cycle from 1d monthly time series (modified from 
;contributed.ncl/rmAnnCycle1D)
;
; The time dimension must be a multiple of 12
; Usage:
;        xAnom = rmAnnCycle1D (x)
;        x     = rmAnnCycle1D (x)


 undef("rmClm1D_smth")
function rmClm1D_smth (x[*]:numeric)
local ntim, nmo, nmos, aveMonth, xAnom, yr
begin
  ntim = dimsizes(x)

  nmos = 12
;  modCheck ("rmAnnCycle1D", ntim, nmos)   ; error check
;
; Compute all 12 monthly averages first. 
;
  aveMonth = new(nmos,typeof(x),getFillValue(x))
  do nmo=0,nmos-1
     aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
  end do

  ;smooth

  aveMonth = runave(aveMonth, 3, -1)


;
; remove the monthly means from each year
;
  xAnom = x          ; variable to variable copy [meta data]
  do yr=0,ntim-1,nmos
     xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
  end do
  xAnom@long_name  = x@long_name + ": Anomalies"
  xAnom@anomaly_op_ncl = "function rmAnnCyc1D "

  return (xAnom)
end

;*********************************
; rmclmDay1D
;***************************************
;remove daily climatology from 1d data

undef("rmClmDay1D")
 function rmClmDay1D(x[*]:numeric, yrst[1]:integer, yren[1]:integer)

;x    - 1d variable
;yrst - first year of data
;yren = final year of data

 local nyrs, DayYr, i, j, cnt, dInd, clm, out, nyrs, cen, cf, xbar

 begin

   yrs   = ispan(yrst, yren, 1)

   ;convert to 2-d array
   nyrs  = dimsizes(yrs)
   DayYr = new((/366, nyrs/), typeof(x))
   
   ;calculate climatology

   ;index of non-leap days
   dInd      = new(365, integer)
   dInd(:58) = ispan(0,58,1)
   dInd(59:) = ispan(60,365,1)

   cnt  = 0   ;counter
   do i = 0, nyrs-1

     if (isleapyear(yrs(i))) then
       
       cen           = cnt+365
       DayYr(:,i)    = x(cnt:cen)
     else

       cen           = cnt+364
       DayYr(dInd,i) = x(cnt:cen)
     end if
     cnt             = cen+1

   end do

   
   clm      = dim_avg(DayYr)

   ;smooth
   cf       = ezfftf(clm)     ;forward fourier transform
   cf(:,2:) = 0.              ;only retain annual and seim-annual harmonics
   xbar     = avg(clm)
   clm      = ezfftb(cf, xbar)
   
   

   ;remove climatology and return to 1d array

   out = x  
   cnt = 0

   do j = 0, nyrs-1


     if (isleapyear(yrs(j))) then

       cen          = cnt+365
       out(cnt:cen) = DayYr(:,j) - clm
     else

       cen          = cnt+364 
       out(cnt:cen) = DayYr(dInd,j) - clm(dInd)
     end if
     
     cnt            = cen+1
     
   end do

   return(out)

 end 




;******************************
;ismonotonic
;*****************************
;Checks a 1-d array to see if it's monotonically increasing
;

 undef("ismonotonic")
 function ismonotonic(x[*]:numeric)

 local n, diff

 begin

   n = dimsizes(x)
   
   diff = x(1:) - x(:n-2)

   if (any(diff.le.0.)) then
     return(False)
   else
     return(True)
   end if

 end


;***********************
;copy_varDimNames
;*******************
;copies dimension names from var1 to var2

 undef("copy_varDimNames")
 procedure copy_varDimNames (var1, var2)

 local dnams, i

 begin

   dnams = getVarDimNames(var1)
   
   do i = 0, dimsizes(dnams)-1     
     var2!i = dnams(i)
   end do

 end


  
;***************************
;CellAverage
;****************************

;Returns gridpoint cosine-weighted average - used for reducing spatial resolution
 undef("CellAverage")
 function CellAverage (x[*][*][*]:numeric, xres[1]:numeric, yres[1]:numeric,lat[*]:numeric)


; compute "block" cosine wgted averages
; missing values [ie: x@_FillValue ] automatically ignored

; x   - input array with dimensions time, lat, lon
; xres - output lon resolution (degrees)
; yres - output lat resolution (degrees)
; lat - 1d array of input lat coords
 
 local dimx, jlat, ilon, dnams, rad, clat, wgt, nlat, mlon, y, nl, jl, ml, il, olat, olon


begin



; input dimension sizes
  dimx = dimsizes(x)      
  jlat = dimx(1)
  ilon = dimx(2)

  ;input dimension names
  dnams = getvardims(x)
 
  rad  = 4.0*atan(1.0)/180.
  clat = cos(lat*rad)

  ;conform lat. weights
  wgt  = new (dimx, typeof(x)) 
  wgt  = conform(wgt, clat, (/1/))


; output dimension sizes
 
  nlat    = floattoint(ceil( (lat(jlat-1) - lat(0)) /yres ) )
  dimx(1) = nlat
  mlon    = floattoint(ceil( (x&$dnams(2)$(ilon-1) - x&$dnams(2)$(0))/xres ))
  dimx(2) = mlon
  y       = new (dimx, typeof(x) )

  nby     = jlat/nlat
  nbx     = ilon/mlon



 ; cosine wgted average: sum(x*wgt)/sum(wgt)
  
  nl = -1
  do jl=0,jlat-1,nby
     nl = nl+1
     ml = -1
    do il=0,ilon-1,nbx

       ml = ml+1
       y(:,nl,ml) = sum(x(:,jl:jl+nby-1,il:il+nbx-1)*wgt(:,jl:jl+nby-1,il:il+nbx-1))/sum(wgt(:,jl:jl+nby-1,il:il+nbx-1))

    end do
  end do

  ;coords
  y!0          = dnams(0)
  y!1          = dnams(1)
  y!2          = dnams(2)

  olat         = fspan(lat(0), lat(jlat-1),nlat)
  olat@units   = lat@units

  olon         =  fspan(x&$dnams(2)$(0), x&$dnams(2)$(ilon-1), mlon) 
  olon@units   = x&$dnams(2)$@units

  y&$dnams(0)$ = x&$dnams(0)$
  y&$dnams(1)$ = olat
  y&$dnams(2)$ = olon

  return(y)

end 
 


;***************************
;CellAverage_nowgt
;****************************

;Returns gridcell average - used for reducing spatial resolution
 undef("CellAverage_nowgt")
 function CellAverage_nowgt (x[*][*][*]:numeric)


; compute grid cell averages
; missing values [ie: x@_FillValue ] automatically ignored

; x   - input array with at least 2 dimensions. Rightmost 2 dims are averaged
 

 local dimx, rank, jlat, ilon, dnams, rad, clat, wgt, nlat, mlon, y, nl, jl, ml, il, olat, olon


begin



; input dimension sizes
  dimx = dimsizes(x)
  rank = dimsizes(dimx)

  if (rank.lt.2.or.rank.gt.4) then
    print("ncl_funcs.ncl - CellAverage_nowgt. x must have 2-4 dimensions")
    exit
  end if
 
  jlat = dimx(rank-2)
  ilon = dimx(rank-1)

  ;input dimension names
  dnams = getvardims(x)



 ; average

  if (rank.eq.2) then
    tmp = 0.5 * (x(1:,:)+x(:ilat-2,:))
    out = 0.5 * (tmp(:,1:) + tmp(:,ilon-2))
  end if

  if (rank.eq.3) then
    tmp = 0.5 * (x(:,1:,:)+x(:,:ilat-2,:))
    out = 0.5 * (tmp(:,:,1:) + tmp(:,:,ilon-2))
  end if

  if (rank.eq.4) then
    tmp = 0.5 * (x(:,:,1:,:)+x(:,:,:ilat-2,:))
    out = 0.5 * (tmp(:,:,:,1:) + tmp(:,:,:,ilon-2))
  end if


  ; dim names
  do i = 0, rank-1
    out!i = dnams(i)
  end do

;coords
  if (rank.gt.2) then

    do i = 0, rank-3
      out&$dnams(i)$ = x&$dnams(i)$
    end do
  end if


  
  lat                 = x&$dnams(rank-2)$
  latAve              = lat(1:)
  latAve              = 0.5*(lat(1:)+lat(:ilat-2) )
  out&$dnams(rank-2)$ = latAve


  lon                 = x&$dnams(rank-1)$
  lonAve              = lon(1:)
  lonAve              = 0.5*(lon(1:)+lat(:ilon-2) )
  out&$dnams(rank-1)$ = lonAve



  return(out)

end 


;********************************
;UT2hourfrac
;****************************

;This function converts a UT time (hhmm) to hourly fraction.
;Useful for time-dependant weight funcs

 undef("UT2hourfrac")
 function UT2hourfrac(t[*]:numeric)

;t - 1-d array of 24hour clock times

 local tInt, hr, mn2dec, out

begin


 if (isfloat(t)) then
   tInt   = floattoint(t)
 else
   tInt   = t
 end if


 ;get hours (ie 0535 returns 5, 1200 returns 12 etc)

   hr = floor(tInt/100.)

 
 ;get mins and convert to fraction of one hour (0536 returns 36/60 = 0.6)

   mn2dec = (tInt%100)/60.


 ;add both

   out = hr + mn2dec

   return(out)

 end

;****************************************
;closest_val_anyord
;***************************************

;taken from Example 2 of ncl closest_val entry

 undef("closest_val_anyord")
function closest_val_anyord(x:numeric, xVal[1]:numeric)

local xDif, xAbsDif, xMinVal, iClose

begin
  xDif    = xVal-x
  xAbsDif = fabs(xDif)
  iClose  = minind(xAbsDif)
  return(iClose)   ; original doc says "first occurence"
end


;***********************************************
;local_max_1d
;********************************************

 undef("local_max_oned")
function local_max_oned(x[*]:numeric, cyclic:logical, del[1]:numeric)

; imx = ind((x(1:nx-2)-x(0:nx-3)).gt.del .and. \
; (x(1:nx-2)-x(2:nx-1)).gt.del) + 1
; if (.not.ismissing(imx)) then
; return(x(imx))
; end if

begin
  if (del.lt.0.0) then
      print("local_max_1d: del must be .ge. 0.0")
      print("local_max_oned: del ="+del)
      exit
  end if

  nx = dimsizes(x)
  imx = new( nx-2, "integer")

  i = -1
  do n=1,nx-2
     if ((x(n)-x(n-1)).ge.del .and. (x(n)-x(n+1)).ge.del) then
         i = i+1
         imx (i) = n
     end if
  end do

  if (i.gt.-1) then
    vals = (ind(.not.ismissing(imx)))
    return(imx(vals))
  else
    return(-999)
  end if


end


;**********************************
;rmv_insuf
;*********************************

;sets data in input array to _FillValue where number of observations is below a specified minimum

  undef("rmv_insuf")
  procedure rmv_insuf(x:numeric, nval:integer, minobs[1]:integer)

;x      - input array of any dimension
;nval   - integer array of number of observations, must be same size as x
;minobs - scalar integer of minimum allowable observation size

  local  nval1d, x1d, nod
  
  begin

	nval1d   = ndtooned(nval)
	nod      = ind(nval1d.lt.minobs)
	x1d      = ndtooned(x)
 
	x1d(nod) = x@_FillValue		
	x        = onedtond(x1d, dimsizes(x))

  end



;**********************************************
;equiv sample size (extends existing func beyond AR(1) processes)
;******************************************

;If this is required see sst_ar_test.ncl.  (von Storch and Zwiers (1995))

  undef("equiv_size")
  function equiv_size (x:numeric,  mxlag[1]:integer, sig[1]:numeric)


;x 	: input array of any dimension. Rightmost dimension should be time. missing values are allowed
;mxlag	: maximum lag to which analysis should be performed. For an AR(p) process, this will be p.
;sig. level : (<1.)


;Return value will be an integer array of size x, less rightmost dimension (time)

   local xDim, xRk, acf, acf1d, mxlag, N, N1d, Nx, nZero, lag, term1, Ne, i, p

begin



;convert to 2d array

  xDim = dimsizes(x)
  xRk  = dimsizes(xDim)  ;rank

	
;autocorrelation func - NB 'acf'includes the (useless) lag-0 autocorr!
   acf  =  abs( esacr(x, mxlag) )    ; dims are (spatial, lag)


;number of non-missing values
   N   = dim_num(.not.ismissing(x))

  if (any(N.lt.30)) then
    print("equiv_size warning: input data has less than 30 data points")
    print("equiv_size is not an appropriate estimator of sample size!")
  end if

  ;add lag dimension to N
   Nx = new(dimsizes(acf), integer)      
   if (xRk.lt.2) then 			;1d case	
   	Nx  = conform(Nx, N, -1)
   else
        Nx  = conform(Nx, N, ispan(0, xRk-2, 1))
   end if	



;correction to acf due to esacr bug (ref Dennis Shea)

   nZero = num(N.eq.0)



   if (nZero.gt.0) then
	
     if (xRk.lt.2) then
       print("ncl_funcs.ncl/equiv_size : input array is empty")
       exit
     else
 
       N1d      = ndtooned(Nx)
       i        = ind(N1d.eq.0)
       acf1d    = ndtooned(acf)
       acf1d(i) = x@_FillValue
       N1d(i)   = Nx@_FillValue   ;this avoids division by zero
       
       acf      = onedtond(acf1d, dimsizes(acf) )
       Nx       = onedtond(N1d, dimsizes(Nx))
     end if        
     
   end if

   
   lag     = conform(acf, ispan(0, mxlag, 1), (/xRk-1/) )
 
   acf = where(lag.eq.0, 0., acf)  ;set lag-zero to 0
	
;remove insig acfs
   if (xRk.eq.1) then
     p = rtest(acf,N,0)
   else
     p   = rtest(acf,conform(acf,N,ispan(0,xRk-2,1)),0)
   end if
   acf = where(p.gt.sig, 0., acf)


;summed term NB! corr. term of 1 is added to accnt for lag-0 acf!
   term1 = 2*dim_sum( (1- (lag/Nx)) * acf )

;check for negative values
   if (any(term1.lt.0)) then
     print("ncl_funcs.ncl/equiv_size : negative values found")
     exit
   end if 	

   	
   Ne   = toint(  N/(1+term1)  )

   return(Ne)

end




;***********************************************
;cfd_1d
;*********************************************

;Performs centered finite differencing on a 1d, non-cyclic array.
;endpoints are set to _FillValue, so that output is same size as input

  undef("cfd_1d")
  function cfd_1d(x[*]:numeric)

  local ndims

  begin

    ndims = dimsizes(dimsizes(x))
    
    if (ndims.gt.1) then
      print("error in cfd_1d - array has >1 dimension")
      exit
    end if

    ndims = dimsizes(x)

    xGrad = new(ndims, typeof(x))
    
    xGrad((/0,ndims-1/)) = x@_FillValue

    do  i = 1, ndims-2

      xGrad(i) = x(i+1) - x(i-1)

    end do

    return(xGrad)

  end


  


;******************************************************************************
;ttest_1samp
;***************************************************************

;Performs one-sample Student's t-test, typically for paired-difference test

;prototype
  undef("ttest_1samp")
  function ttest_1samp(x:numeric, var:numeric, N:integer, c[1]:numeric, tval[1]:logical)

;description

;x	: array of any dimension of means
;var	: array of variances, must be same dimensions as x   	
;N      : scalar or array of sam size as x of sample size
;c      : scalar value against which x is compared. Typically 0.
;tval	: if True, t values are returned, if false probabilities are returned. 

;Null hypothesis is that x - c = 0

  local xdim, vdim, ndim, t, df, prob, uplim, beta_a, beta_b, ii, var1d

  begin



	;test that x, var are same dimensions

	xdim = dimsizes(x)
	vdim = dimsizes(var)

	if (any(xdim.ne.vdim)) then
	   print("fatal: ttest_1samp. Mean and variance arrays must have same dimensions")
	   exit
	end if		

	;test that N is scalar or correctly sized array

	ndim = dimsizes(N)
	

	if (sum(ndim).gt.1.and.any(ndim.ne.vdim)) then
		print("fatal: ttest_1samp. Sample size must be a scalar, or have same dimensions as x")
		exit
	end if

        if (any(N.le.2)) then
		print("fatal: ttest_1samp. Sample size must be greater than 2")
		exit          
        end if
		 


;calculate t-value


        if (any(var.le.0)) then

          var1d     = ndtooned(var)
          ii        = ind(var1d.le.0)
          var1d(ii) = var@_FillValue
          var       = onedtond(var1d, vdim)
        end if

	t     = (x - c) / (( var/N)^0.5 )



        if (tval.eq.True) then

          return(t)
                
        else
    
          if (num(N).eq.1) then
            df     = conform(t, N-2, -1) 
          else
            df     = N-2
          end if
            uplim  = df/(df+t^2)
            beta_a = df/2.
            beta_b = conform(t, 0.5, -1)
            prob   = betainc(uplim , beta_a, beta_b)

            return(prob) 		
          end if
	

  end

