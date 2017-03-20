load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01011993_31121997_T.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01021998_31031998.nc","r")
	data = addfile("data/dataset-global-reanalysis-phy-001-025-31121992_31122003_T&U&V.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-30121992_29122015_region2.nc","r")
	
	; temperature time series
	TAB = str_get_tab()
		

; *********************** DATA *****************************************************************

	;sws 	= short2flt(data->salinity(:,0,:,:))	; sea water salinity
	temp 	= short2flt(data->temperature(:,0,:,:)) ; temperature
	;v 	= short2flt(data->v(:,0,:,:))		; northward velocity
	;u 	= short2flt(data->u(:,0,:,:))		; eastward velocity
	lat 	= data->latitude			; latitude
	lon 	= data->longitude			; longitude		
	time 	= cd_calendar((data->time),0)		; time
	;dep 	= data->depth				; depth
	;spd 	= sqrt(u^2+v^2)				; scalar velocity

	temp 		= temp - 273.15	; K->C
	temp@units 	= "(C)"
	temp@unit_long 	= "Celsius"
	temp@lat 	= lat
	temp@lon 	= lon

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

; *********************** PLOT TEMPERATURE VS TIME *******************************************
	wks 					= gsn_open_wks("x11","Tvst")

	res 					= True

	;wks 					= gsn_open_wks("x11",file_title)
	res@gsnRightString			= ""+lat_UpB+"~S~o~N~N to "+lat_LoB+"~S~o~N~N, "+lon_LeB+"~S~o~N~E to "+lon_RiB+"~S~o~N~E"
	res@gsnRightStringFontHeightF 		= 0.015
	res@gsnRightStringOrthogonalPosF 	= 0.007
	res@gsnLeftStringFontHeightF 		= 0.015
	res@gsnLeftStringOrthogonalPosF 	= 0.007
	;res@sfXArray 				= lon
	;res@sfYArray 				= lat
	res@tmXBLabelFontHeightF 		= 0.014
	res@tmYLLabelFontHeightF 		= 0.014

	ndays 					= dimsizes(dd)
	mean_T 					= new(ndays,float)
	dateF 					= new(ndays,float)

	; vector date x-axis
	dateF(0)=1
	do n=1,ndays-1
		dateF(n) = dateF(n-1)+1
	end do

	res@tiMainString 			= "Mean Sea Surface Temperature"
	res@tiXAxisString 			= ""
	res@tiYAxisString 			= "Temperature (~S~o~N~C)"
	res@tiYAxisFontHeightF 			= 0.018
	res@tiXAxisString 			= "Days since "+dd(0)+"-"+mm(0)+"-"+yyyy(0)+""
	res@tiMainOffsetYF 			= 0.03
	res@xyMarkLineModes 			= "Markers"
	res@xyMarkerColor 			= "black"
	res@xyMarkers 				= 16
	res@xyMarkerSizeF 			= 0.005

	;res@tmXBMode        = "Manual"	
	;res@tmXBTickStartF  = yyyy(0)
	;res@tmXBTickEndF    = yyyy(dimsizes(yyyy)-1)
	;res@tmXBTickSpacingF= 1

	mean_T 					= dim_avg_n(dim_avg_n(temp(:,LoB:UpB,LeB:RiB),1),1)

	res@gsnLeftString  			= " From "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ range=["+min(mean_T)+","+max(mean_T)+"]"

	
	plot = gsn_csm_xy(wks,dateF,mean_T,res)



end
