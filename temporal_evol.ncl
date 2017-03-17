load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

	data = addfile("data/dataset-global-reanalysis-phy-001-025-01011993_31121997_T.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phy-001-025-01021998_31031998.nc","r")
	
	; temperature time series
	
	

; *********************** DATA *************************************************************

	;sws = short2flt(data->salinity(:,0,:,:))	; sea water salinity
	temp = short2flt(data->temperature(:,0,:,:)) 	; temperature
	;v = short2flt(data->v(:,0,:,:))		; northward velocity
	;u = short2flt(data->u(:,0,:,:))		; eastward velocity
	lat = data->latitude				; latitude
	lon = data->longitude				; longitude		
	time = cd_calendar((data->time),0)		; time
	dep = data->depth				; depth
	;spd = sqrt(u^2+v^2)				; scalar velocity

	temp = temp - 273.15	; K->C
	temp@units = "(C)"
	temp@unit_long = "Celsius"
	temp@lat = lat
	temp@lon = lon

	yyyy	=floattointeger(time(:,0))		; full year
	mm	=floattointeger(time(:,1))		; month
	dd	=floattointeger(time(:,2))		; day of month


; *********************** INITIALIZATION *************************************************************

		

	res 					= True

	;wks 					= gsn_open_wks("x11",file_title)

	res@cnMonoLineColor 			= False
	res@cnFillOn 				= True			; Fill with color
	res@cnLinesOn 				= False			; No contour lines

	res@pmTickMarkDisplayMode 		= "Always" 	; Default lat/lon labels
	;res@gsnMaximize 			= "True"		; Maximize plot in frame
	res@gsnAddCyclic          		= False       ; data already has cyclic point

	res@tiMainOffsetYF 			= -0.03
		
	res@gsnLeftStringFontHeightF 		= 0.015
	res@gsnLeftStringOrthogonalPosF 	= 0.002
	res@gsnRightStringFontHeightF 		= 0.015
	res@gsnRightStringOrthogonalPosF 	= 0.006
  	;res@gsnScalarContour 			= True	; contours desired
		
	; Zoom in on map area that we are interested in
	res@mpLimitMode      			= "Corners"
  	res@mpLeftCornerLatF 			= lat(0)
  	res@mpLeftCornerLonF 			= lon(0)
  	res@mpRightCornerLatF 			= lat(dimsizes(lat)-1)
  	res@mpRightCornerLonF 			= lon(dimsizes(lon)-1)
	res@mpOutlineBoundarySets		= "GeophysicalAndUSStates"
	res@mpDataBaseVersion 			= "HighRes"

	res@lbOrientation        		= "Vertical"
	res@lbTitlePosition      		= "Right"                           ; title location
  	res@lbTitleDirection     		= "Across"                          ; letter angle
  	res@lbTitleAngleF        		= 90.                               ; title angle
  	res@lbTitleFontHeightF   		= 0.02                              ; font height
	res@lbTitleOffsetF	 		= 0.3

	res@sfXArray 				= lon
	res@sfYArray 				= lat
	res@pmLabelBarDisplayMode 		= "Always"
	res@lbPerimOn             		= False
	res@cnInfoLabelOn 			= False
	res@tmXBLabelFontHeightF 		= 0.018
	res@tmYLLabelFontHeightF 		= 0.018
	
	ndays 					= dimsizes(dd)
	mean_T 					= new(ndays,"float")
	dateF 					= new(ndays,float)

	; vector date x-axis
	dateF(0)=1
	do n=1,ndays-1
		dateF(n) = dateF(n-1)+1
	end do

	; resources for time serie graphic
	resE 					= True

	resE@tiMainString 			= "Mean Surface Temperature"
	resE@tiXAxisString 			= ""
	resE@tiYAxisString 			= "Temperature (~S~o~N~C)"
	resE@tiXAxisString 			= "Days since "+dd(0)+"-"+mm(0)+"-"+yyyy(0)+""
	resE@xyMarkLineModes 			= "Markers"
	resE@xyMarkerColor 			= "black"
	resE@xyMarkers 				= 16
	resE@xyMarkerSizeF 			= 0.01

	resE@trYMinF     			= 26.5		; set minimum Y-axis value
  	resE@trYMaxF     			= 28		; set maximum Y-axis value
  	;resE@trXMinF     			= min(years)-1        ; set minimum X-axis value
  	;resE@trXMaxF     			= max(years)+1        ; set maximum X-axis value

;*********************** LOOP through every day to create a different file ***************************
		
	
	system("mkdir files_for_animation")
	path = "/home/surf/Desktop/scripts_elise/files_for_animation/"

	n=0
	do n=0,ndays-1
		plot = new(2,graphic)
		wks = gsn_open_wks("png",path+"animate"+sprinti("%03i",n)) ; animate_00.pdf

	; *********************** PLOT TEMPERATURE ***************************************************

		;wks 				= gsn_open_wks("x11","T")
		resT 				= res		
			
		cmap 				= read_colormap_file("BkBlAqGrYeOrReViWh200")
		resT@cnFillPalette 		= cmap(18:177,:) ;ncolor=200

		;resT@tiMainString 		= "Sea Surface Temperature (C)"
		resT@gsnLeftString  		= ""+mm(n)+"/"+dd(n)+"/"+yyyy(n)+" ~C~ z="+dep(0)+"m ~C~ range=["+min(temp(n,:,:))+","+max(temp(n,:,:))+"]"
  		resT@gsnRightString  		= "~S~o~N~C"; right title
		resT@lbTitleString		= "Temperature (~S~o~N~C)" ; bar title
		
		;Select our own contour levels.
  		resT@cnLevelSelectionMode 	= "ManualLevels"
  		resT@cnLevelSpacingF      	= 0.2
		resT@cnMinLevelValF 		= 14.0
		resT@cnMaxLevelValF 		= 32.0
			
		plot(0) 			= gsn_csm_contour_map(wks,temp(n,:,:),resT)

	; *********************** PLOT TEMP VS TIME **************************************************

		;wks 				= gsn_open_wks("x11","temperature evolution")
		mean_T(n) 			= avg(temp(n,:,:))
		plot(1) 			= gsn_csm_xy(wks,dateF,mean_T,resE)

		gsn_panel(wks,plot,(/2,1/),False) ; create panel

		
	end do

	; convert images to animated GIF
	system("convert -delay 10 files_for_animation/animate*.000003.png animations/animated_T.gif")
	system("rm -rf files_for_animation/") ; delete folder


end
