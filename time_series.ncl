load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01011993_31121997_T.nc","r")
	data = addfile("data/dataset-global-reanalysis-phy-001-025-30121992_29122015_-5.-12N_8.12E.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-31121992_31122003_T&U&V.nc","r")
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


	lat_UpB =  -8				; coordinates of the analyzed section
	lat_LoB =  -12				; in degrees east and north
	lon_LeB =   8
	lon_RiB =   12
	UpB = ind(lat.eq.lat_UpB)		; Upper Boundary
	LoB = ind(lat.eq.lat_LoB)		; Lower Boundary
	LeB = ind(lon.eq.lon_LeB)		; Left  Boundary
	RiB = ind(lon.eq.lon_RiB)		; Right Boundary

; *********************** PLOT TEMPERATURE VS TIME *******************************************
	file_title 				= ""+dd(0)+"-"+mm(0)+"-"+yyyy(0)+"_"+dd(dimsizes(dd)-1)+"-"+mm(dimsizes(dd)-1)+"-"+yyyy(dimsizes(dd)-1)+"_"+lat_UpB+"N."+lat_LoB+"N."+lon_LeB+"E."+lon_RiB+"E"
	
	wks 					= gsn_open_wks("ps","msst_"+file_title)

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
	res@tiMainFontHeightF 			= 0.018
	res@tiXAxisString 			= ""
	res@tiYAxisString 			= "Temperature (~S~o~N~C)"
	res@tiYAxisFontHeightF 			= 0.018
	;res@tiXAxisString 			= "Days since "+dd(0)+"-"+mm(0)+"-"+yyyy(0)+""
	res@tiMainOffsetYF 			= 0.02
	res@xyMarkLineModes 			= "Markers"
	res@xyMarkerColor 			= "black"
	res@xyMarkers 				= 16
	res@xyMarkerSizeF 			= 0.005


; --------------- TIME TICK FUNCTION -------------------------------------------
	maxlen = 40
	posx = new(maxlen,integer)
	i=0
	j=0
	mon = (/"jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"/)
	if dimsizes(dd).ge.(365*3) then ; year tick
		tickx = new(maxlen,integer)
		do while (i.le.dimsizes(dd)-1)
			if dd(i).eq.1 .and. mm(i).eq.1 .or. i.eq.(dimsizes(dd)-1)
				posx(j) = i
				tickx(j) = yyyy(i)
				if i.eq.(dimsizes(dd)-1) then tickx(j) = yyyy(i)+1 end if
				j = j+1
			end if
			i=i+1
		end do
	else 		; month-year tick
		tickx = new(maxlen,string)
		do while (i.le.dimsizes(dd)-1)
			if dd(i).eq.1 .or. i.eq.(dimsizes(dd)-1)
				posx(j) = i
				y_tick = yyyy(i) - floattoint(yyyy(i)/100)*100
				if y_tick.le.10 then zero="0"
				else zero=""
				end if
				tickx(j) = ""+mon(mm(i)-1)+""+zero+""+y_tick+""
				if i.eq.(dimsizes(dd)-1) then
					k=0
					if mm(i).eq.12 then
						k=-11
						y_tick = y_tick+1
					end if
					tickx(j) = ""+mon(mm(i)-1+k)+""+zero+""+y_tick+""
				end if
				j = j+1
			end if
			i=i+1
		end do
	end if

	res@tmXBMode             		= "Explicit"	; Tick Marks/Labels
	res@tmXBValues           		= posx(0:j)
	res@tmXBLabels           		= tickx(0:j)
	
	; ---------------- END TIME TICK FUNCTION ---------------------------------------

	;res@trXMinF     = 0        ; set minimum X-axis value
	res@trXMaxF     = dimsizes(dd)-1        ; set maximum X-axis value
	res@gsnMaximize 			= True
	res@vpHeightF 				= 0.3        ; Changes the aspect ratio
	res@vpWidthF  				= 1.3

	mean_T 					= dim_avg_n(dim_avg_n(temp(:,LoB:UpB,LeB:RiB),1),1)

	res@gsnLeftString  			= " From "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ range=["+min(mean_T)+","+max(mean_T)+"]"

	
	plot = gsn_csm_xy(wks,dateF,mean_T,res)



end
