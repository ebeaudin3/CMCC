load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin

	data = addfile("data/dataset-global-reanalysis-phy-001-025-01091994_31081996.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01021998_31031998.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-31121992_31122003_T&U&V.nc","r")
	
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
	
	file_title 				= ""+dd(0)+"-"+mm(0)+"-"+yyyy(0)+"_"+dd(dimsizes(dd)-1)+"-"+mm(dimsizes(dd)-1)+"-"+yyyy(dimsizes(dd)-1)+""
	
	wks 					= gsn_open_wks("x11","vec_"+file_title)

	mean_u 					= transpose(dim_avg_n(dim_avg_n(u(:,LoB:UpB,LeB:RiB),1),1))
	mean_v 					= transpose(dim_avg_n(dim_avg_n(v(:,LoB:UpB,LeB:RiB),1),1))
	mean_spd				= transpose(dim_avg_n(dim_avg_n(spd(:,LoB:UpB,LeB:RiB),1),1))


	; ------- temporal average ------- ;
	
	time_average = 3 			; 0: daily, 1: monthly, 2: yearly, 3: bi-mensual, 4: 4-annual
	row = 3
	col = dimsizes(mean_u(:,0))
	U  = new((/row,col/),float)
	V  = U
	stick_pos = 0

	if (time_average.ge.1) then
		if (time_average.eq.1) then 		; monthly
			val_ind = ind(dd.eq.1)
		else if (time_average.eq.2) then	; yearly
			val_ind = ind(dd.eq.1 .and. mm.eq.1)
		else if (time_average.eq.3) then	; bi-mensual
			val_ind = ind((dd.eq.1 .or. dd.eq.15))
		else if (time_average.eq.4) then	; 4-annual
			val_ind = ind(dd.eq.1 .and. (mm.eq.1 .or. mm.eq.4 .or. mm.eq.7 .or. mm.eq.10))
		end if
		end if
		end if
		end if
		i=0
		do while i.le.(dimsizes(val_ind)-2)
			U(stick_pos,val_ind(i)) = avg(mean_u(val_ind(i):val_ind(i+1),0))
			V(stick_pos,val_ind(i)) = avg(mean_v(val_ind(i):val_ind(i+1),0))
			i=i+1
		end do
	else if (time_average.eq.0) then		; daily
		U(stick_pos,:)	= mean_u(:,0)
		V(stick_pos,:) 	= mean_v(:,0)

	end if
	end if

; ---------------------- graphic resources ------------------------------------------ ;
	resV 					= True

	resV@vpXF 				= 0.04
	resV@vpWidthF 				= 0.94
	resV@vpHeightF 				= (2./3.)*resV@vpWidthF
	
	;resV@tiXAxisString 			= "Days since "+dd(0)+"-"+mm(0)+"-"+yyyy(0)+""
	resV@tiMainString 			= "Mean Current Magnitude and Direction ~C~"+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+"~C~"+lat_UpB+"~S~o~N~N to "+lat_LoB+"~S~o~N~N, "+lon_LeB+"~S~o~N~E to "+lon_RiB+"~S~o~N~E "
	resV@tiMainFontHeightF 			= 0.015
	resV@tiXAxisFontHeightF 		= 0.15
	resV@tiMainOffsetXF			= -0.2
	resV@tiMainOffsetYF			= 0.04

	resV@vcRefMagnitudeF 			= 0.2		; arrows resources
	resV@vcRefLengthF 			= 0.08
	resV@vcGlyphStyle 			= "LineArrow"
	resV@vcMapDirection 			= False
	resV@vcLineArrowThicknessF 		= 0.9
	resV@vcPositionMode 			= "ArrowTail"
	resV@vcRefAnnoOn 			= True		; Turns on reference vector box (default = True)
	resV@vcRefAnnoOrthogonalPosF 		= .15		; Moves reference vector box inside the plot
	resV@vcRefAnnoPerimOn 			= False
	resV@vcRefAnnoString1 			= ""+resV@vcRefMagnitudeF+" m/s"
	resV@vcRefAnnoString2			= ""
	;resV@vcRefAnnoSide			= "Up"

	resV@tmYLMode             		= "Explicit"	; Tick Marks/Labels
	resV@tmYLValues           		= (/0,1,2/)
	resV@tmYLLabels           		= (/"","",""/)
	resV@tmYLLabelAngleF      		= 90.
	resV@tmYLLabelJust        		= "CenterCenter"
	resV@tmYLLabelFontHeightF 		= 0.015
	resV@tmXBLabelFontHeightF		= 0.015

	resV@tmXBBorderOn      = True		; background
	resV@tmXBOn            = True
	resV@tmXTBorderOn      = False
	resV@tmXTOn            = False
	resV@tmYLBorderOn      = False
	resV@tmYLOn            = False
	resV@tmYRBorderOn      = False
	resV@tmYROn            = False


	; --------------- TIME TICK FUNCTION -------------------------------------------
	maxlen = 40
	posx = new(maxlen,integer)
	i=0
	j=0
	mon = (/"jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"/)
	if dimsizes(dd).ge.(365*2) then ; year tick
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

	resV@tmXBMode             		= "Explicit"	; Tick Marks/Labels
	resV@tmXBValues           		= posx(0:j)
	resV@tmXBLabels           		= tickx(0:j)
	
	; ---------------- END TIME TICK FUNCTION ---------------------------------------

	
	; these four resources allow the user to stretch the plot size, and
	; decide exactly where on the page to draw it.
	 resV@gsnMaximize = True
	 resV@vpHeightF = 0.1        ; Changes the aspect ratio ;0.6-2.5
	 resV@vpWidthF  = 0.8		; bien pour 10 ans: 0.1-0.8

	;restick = True
	;restick@ttmFormat = "%d%c"
	;time_axis_labels(time_d,resV,restick) ; call the formatting procedure

	vector = gsn_vector(wks,U,V,resV)


end
