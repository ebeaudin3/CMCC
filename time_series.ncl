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
		

; ***************************** DATA **********************************

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


	lat_UpB =  -5				; coordinates of the analyzed section
	lat_LoB =  -10				; in degrees east and north
	lon_LeB =   10
	lon_RiB =   12
	UpB = ind(lat.eq.lat_UpB)		; Upper Boundary
	LoB = ind(lat.eq.lat_LoB)		; Lower Boundary
	LeB = ind(lon.eq.lon_LeB)		; Left  Boundary
	RiB = ind(lon.eq.lon_RiB)		; Right Boundary


	;--------- El Nino indices ----------;
	data_nino = asciiread("../noaa/sstoi.indices",-1,"float")
	nrow = 422
	elNino = new((/nrow,10/),"float")
	njump = 5
	i=0
	do while i.le.(nrow-1)
		ind_s = i*10+njump
		ind_e = i*10+9+njump

		if ind_e.gt.dimsizes(data_nino) then
			ind_e = dimsizes(data_nino)
		end if
		elNino(i,:) = data_nino((ind_s):(ind_e))
		i=i+1	
	end do
	
	elNino(:,0) = floattoint(elNino(:,0)) ; years
	elNino(:,1) = floattoint(elNino(:,1)) ; months

	; display elNino 3.4 data on same length vector as RA data
	Nino = new((/dimsizes(dd),10/),"float")
	i=0
	ind_2=0
	do while i.le.(dimsizes(dd)-1)
		ind_1 = ind(yyyy(i).eq.elNino(:,0) .and. mm(i).eq.elNino(:,1))
		if ind_1.gt.ind_2 then
			Nino(i+14,:) = elNino(ind_1,:) ; 8: Nino3.4 temperature
			ind_2=ind_1
		end if
		i=i+1
	end do
	;-------- END Nino indices ----------;
	


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

	ndays 					= dimsizes(dd)
	mean_T 					= new(ndays,float)
	dateF 					= new(ndays,float)

	; vector date x-axis
	dateF(0)=1
	do n=1,ndays-1
		dateF(n) = dateF(n-1)+1
	end do

	mean_T 					= dim_avg_n(dim_avg_n(temp(:,LoB:UpB,LeB:RiB),1),1)

; **************************** GRAPHIC RESOURCES ******************************
	file_title 				= ""+dd(0)+"-"+mm(0)+"-"+yyyy(0)+"_"+dd(dimsizes(dd)-1)+"-"+mm(dimsizes(dd)-1)+"-"+yyyy(dimsizes(dd)-1)+"_"+lat_UpB+"N."+lat_LoB+"N."+lon_LeB+"E."+lon_RiB+"E"
	wks 					= gsn_open_wks("png","msst_"+file_title+"Nino3.4_Nino1+2")

	res 					= True

	res@gsnRightString			= ""+lat_UpB+"~S~o~N~N to "+lat_LoB+"~S~o~N~N, "+lon_LeB+"~S~o~N~E to "+lon_RiB+"~S~o~N~E"
	res@gsnRightStringFontHeightF 		= 0.015
	res@gsnRightStringOrthogonalPosF 	= 0.007
	res@gsnLeftStringFontHeightF 		= 0.015
	res@gsnLeftStringOrthogonalPosF 	= 0.007
	;res@sfXArray 				= lon
	;res@sfYArray 				= lat
	res@tmXBLabelFontHeightF 		= 0.014
	res@tmYLLabelFontHeightF 		= 0.014

	;res@tiMainString 			= "Mean Sea Surface Temperature"
	;res@tiMainFontHeightF 			= 0.018
	res@tiXAxisString 			= ""
	res@tiYAxisString 			= "Temperature (~S~o~N~C)"
	res@tiYAxisFontHeightF 			= 0.018
	res@tiMainOffsetYF 			= 0.02


	res@pmLegendDisplayMode    		= "Always"            	; turn on legend
	res@pmLegendSide           		= "Right"              	; Change location of 
	res@pmLegendParallelPosF   		= .95                  	; move units right
	res@pmLegendOrthogonalPosF 		= -0.2                	; move units down
	res@pmLegendWidthF         		= 0.15                	; Change width and
	res@pmLegendHeightF        		= 0.05                	; height of legend.
	res@lgPerimOn              		= False               	; turn off box around
	res@lgLabelFontHeightF     		= .01                 	; label font height
	res@xyExplicitLegendLabels 		= "SST Angola"         	; create explicit labels
	
	
	res@tmXBMode             		= "Explicit"		; Tick Marks/Labels
	res@tmXBValues           		= posx(0:j)
	res@tmXBLabels           		= tickx(0:j)
	
	res@trXMaxF    				= dimsizes(dd)-1        ; set maximum X-axis value
	res@trYMaxF				= 32
	res@gsnMaximize 			= True
	res@vpHeightF 				= 0.3        		; Changes the aspect ratio
	res@vpWidthF  				= 1.3
	
	res@gsnLeftString  			= " From "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ range=["+min(mean_T)+","+max(mean_T)+"]"

	res@xyLineThicknesses 			= 1.5
 	res@xyLineColors      			= "blue"
	

; ------------ overlay el nino 3.4 index ------------------------------
	plotRA = gsn_csm_xy(wks,dateF,mean_T,res)
	
	res@xyExplicitLegendLabels 		= "SST El Nino 3.4"         	; create explicit labels
	res@pmLegendParallelPosF   		= res@pmLegendParallelPosF-0.05; move units right
	res@xyLineThicknesses 			= 1.0
	res@xyLineColors      			= "red"
	
	index_no_msg = ind(.not.ismissing(Nino(:,8)))
	plotNino34 = gsn_csm_xy(wks,dateF(index_no_msg),Nino(index_no_msg,8),res)
	

	res@xyExplicitLegendLabels 		= "SST El Nino 1+2"         	; create explicit labels
	res@pmLegendParallelPosF   		= res@pmLegendParallelPosF-0.05; move units right
	res@xyLineThicknesses 			= 1.0
	res@xyLineColors      			= "green"
	
	plotNino12 = gsn_csm_xy(wks,dateF(index_no_msg),Nino(index_no_msg,4),res)
	

	overlay(plotRA,plotNino34)
	overlay(plotRA,plotNino12)
	draw(plotRA)
	frame(wks)


; ------------------------- Distance between peaks -----------------------------------





end
