load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin

	data = addfile("data/dataset-global-reanalysis-phy-001-025-01012004_31122004.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01021998_31031998.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-31121992_31122003_T&U&V.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-30121992_29122015_region2.nc","r")
	
	; temperature time series
	TAB = str_get_tab()
		

; *********************** DATA *****************************************************************

	;sws 	= short2flt(data->salinity(:,0,:,:))	; sea water salinity
	;temp 	= short2flt(data->temperature(:,0,:,:)) ; temperature
	v 	= short2flt(data->v(:,0,:,:))		; northward velocity
	u 	= short2flt(data->u(:,0,:,:))		; eastward velocity
	lat 	= data->latitude			; latitude
	lon 	= data->longitude			; longitude
	time_d  = data->time		
	time 	= cd_calendar((data->time),0)		; time
	;dep 	= data->depth				; depth
	spd 	= sqrt(u^2+v^2)				; scalar velocity

	;temp 		= temp - 273.15	; K->C
	;temp@units 	= "(C)"
	;temp@unit_long 	= "Celsius"
	;temp@lat 	= lat
	;temp@lon 	= lon

	yyyy	= floattointeger(time(:,0))		; full year
	mm	= floattointeger(time(:,1))		; month
	dd	= floattointeger(time(:,2))		; day of month


	lat_UpB =  2				; coordinates of the analyzed section
	lat_LoB =  -2				; in east and north
	lon_LeB =   2
	lon_RiB =   10
	UpB = ind(lat.eq.lat_UpB)		; Upper Boundary
	LoB = ind(lat.eq.lat_LoB)		; Lower Boundary
	LeB = ind(lon.eq.lon_LeB)		; Left  Boundary
	RiB = ind(lon.eq.lon_RiB)		; Right Boundary


; *********************** PLOT MEAN VELOCITY VS TIME ******************************************
	wks 					= gsn_open_wks("x11","SPDvst")

	mean_u 					= transpose(dim_avg_n(dim_avg_n(u(:,LoB:UpB,LeB:RiB),1),1))
	mean_v 					= transpose(dim_avg_n(dim_avg_n(v(:,LoB:UpB,LeB:RiB),1),1))
	mean_spd				= transpose(dim_avg_n(dim_avg_n(spd(:,LoB:UpB,LeB:RiB),1),1))

	row 					= 3
	col 					= dimsizes(mean_u(:,0))

	U   					= new((/row,col/),float)
	V   					= U
	stick_pos  				= 1
	U(stick_pos,:) 				= mean_u(:,0)
	V(stick_pos,:) 				= mean_v(:,0)
	  
	resV 					= True

	resV@vpXF 				= 0.04
	resV@vpWidthF 				= 0.94
	resV@vpHeightF 				= (2./3.)*resV@vpWidthF
	
	resV@tiXAxisString 			= "Days since "+dd(0)+"-"+mm(0)+"-"+yyyy(0)+""
	resV@tiMainString 			= "Mean Current Magnitude and Direction ~C~"+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+"~C~"+lat_UpB+"~S~o~N~N to "+lat_LoB+"~S~o~N~N, "+lon_LeB+"~S~o~N~E to "+lon_RiB+"~S~o~N~E "
	resV@tiXAxisFontHeightF 		= 0.015
	resV@tiMainFontHeightF 			= 0.018
	resV@tiMainOffsetXF			= -0.2
	resV@tiMainOffsetYF			= 0.04
	resV@tmXBLabelFontHeightF		= 0.015


	resV@vcRefMagnitudeF 			= 0.5
	resV@vcRefLengthF 			= 0.08
	;resV@vcMinFracLengthF 			= 0
	resV@vcGlyphStyle 			= "LineArrow"
	resV@vcMapDirection 			= False
	resV@vcLineArrowThicknessF 		= 2.5
	resV@vcPositionMode 			= "ArrowTail"
	resV@vcRefAnnoOn 			= True		; Turns on reference vector box (default = True)
	resV@vcRefAnnoOrthogonalPosF 		= -.15		; Moves reference vector box inside the plot
	resV@vcRefAnnoPerimOn 			= False
	resV@vcRefAnnoString2 			= ""
	resV@vcRefAnnoSide			= "Right"

	resV@tmYLMode             		= "Explicit"	; Tick Marks/Labels
	resV@tmYLValues           		= (/0,1,2/)
	resV@tmYLLabels           		= (/"","",""/)
	resV@tmYLLabelAngleF      		= 90.
	resV@tmYLLabelJust        		= "CenterCenter"
	resV@tmYLLabelFontHeightF 		= 0.015
	;resV@tmXBMode             		= "Explicit"	; Tick Marks/Labels
	;resV@tmXBValues           		= (/0,150,300/)
	;resV@tmXBLabels           		= (/"","",""/)
	
	resV@sfXArray 	 = time_d
	
	; these four resources allow the user to stretch the plot size, and
	; decide exactly where on the page to draw it.

	 resV@gsnMaximize = True
	 resV@vpHeightF = 0.3        ; Changes the aspect ratio
	 resV@vpWidthF  = 0.9

	;restick = True
	;restick@ttmFormat = "%d%c"
	;time_axis_labels(time_d,resV,restick) ; call the formatting procedure

	vector = gsn_vector(wks,U,V,resV)


end
