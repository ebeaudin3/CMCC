load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

	data = addfile("data/dataset-global-reanalysis-phy-001-025-01012004_31122004.nc","r")
	
	; temperature, salinity, velocity
	
	i=1
	if (i.eq.0) then
		print(data)
	end if

	if (i.eq.1) then

; ----------- DATA -------------------------------------------------

		sws = short2flt(data->salinity); sea water salinity
		temp = short2flt(data->temperature) ; temperature
		v = short2flt(data->v)		; northward velocity
		u = short2flt(data->u)		; eastward velocity
		lat = data->latitude			; latitude
		lon = data->longitude			; longitude		
		time = cd_calendar((data->time),0)	; time
		dep = data->depth			; depth
		spd = sqrt(u^2+v^2)			; scalar velocity

		temp = temp - 273.15	; K->C
		temp@units = "(C)"
		temp@lat = lat
		temp@lon = lon

		yyyy        =floattointeger(time(:,0))	; full year
		mm          =floattointeger(time(:,1))	; month
		dd          =floattointeger(time(:,2))	; day of month

		

; ----------- INITIALIZATION -------------------------------------

		res 				= True

		plot 				= new(3,graphic)
		file_title 			= ""+dd(0)+"-"+mm(0)+"-"+yyyy(0)+"_"+dd(dimsizes(dd)-1)+"-"+mm(dimsizes(dd)-1)+"-"+yyyy(dimsizes(dd)-1)+""
		wks = gsn_open_wks("ps",file_title)

		res@gsnDraw  			= False			; don't draw
		res@gsnFrame 			= False			; don't advance frame		
		res@gsnLeftStringFontHeightF 	= 0.015
		res@gsnLeftStringOrthogonalPosF = 0.002
		res@gsnRightStringFontHeightF 	= 0.015
		res@gsnRightStringOrthogonalPosF= 0.006
  		res@gsnScalarContour 		= True	; contours desired
		;res@gsnMaximize = "True"		; Maximize plot in frame
		res@gsnAddCyclic          	= False       ; data already has cyclic point

		res@cnMonoLineColor 		= False
		res@cnFillOn 			= True			; Fill with color
		res@cnLinesOn 			= False			; No contour lines
		res@cnInfoLabelOn 		= False
		
		res@pmTickMarkDisplayMode 	= "Always" 	; Default lat/lon labels
		res@pmLabelBarDisplayMode 	= "Always"

		res@tiMainOffsetYF 		= -0.03

		; Zoom in on map area that we are interested in
		res@mpLimitMode      		= "Corners"
  		res@mpLeftCornerLatF 		= lat(0)
  		res@mpLeftCornerLonF 		= lon(0)
  		res@mpRightCornerLatF 		= lat(dimsizes(lat)-1)
  		res@mpRightCornerLonF 		= lon(dimsizes(lon)-1)
		res@mpOutlineBoundarySets 	= "GeophysicalAndUSStates"
		res@mpDataBaseVersion 		= "HighRes"

		res@lbOrientation        	= "Vertical"
		res@lbTitlePosition      	= "Right"                           ; title location
  		res@lbTitleDirection     	= "Across"                          ; letter angle
  		res@lbTitleAngleF       	= 90.                               ; title angle
  		res@lbTitleFontHeightF   	= 0.02                              ; font height
		res@lbTitleOffsetF	 	= 0.3
		res@lbPerimOn            	= False
		
		res@sfXArray 			= lon
		res@sfYArray 			= lat

		res@tmXBLabelFontHeightF 	= 0.018
		res@tmYLLabelFontHeightF 	= 0.018

; ----------- PLOT TEMPERATURE -------------------------------------

		;wks = gsn_open_wks ("x11","T")

		resT 				= res		
		
		cmap 				= read_colormap_file("BkBlAqGrYeOrReViWh200")
		resT@cnFillPalette 		= cmap(18:177,:) ;ncolor=200

		;resT@tiMainString 		= "Sea Surface Temperature (C)"
		resT@gsnLeftString  		= " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ z="+dep(0)+"m ~C~ range=["+min(dim_avg_n(temp(:,0,:,:),0))+","+max(dim_avg_n(temp(:,0,:,:),0))+"]"
  		resT@gsnRightString  		= "~S~o~N~C"; right title
		resT@lbTitleString 		= "Temperature (~S~o~N~C)" ; bar title
		
		;Select our own contour levels.
  		resT@cnLevelSelectionMode 	= "ManualLevels"
  		resT@cnLevelSpacingF      	= 0.4
		resT@cnMinLevelValF 		= 12.0
		resT@cnMaxLevelValF 		= 32.4
		;resT@lbLabelStrings 		= (/"12","","","","","14","","","","","16","","","","","18","","","",""\
		;				,"20","","","","","22","","","","","24","","","","","26","","","",""\
		;				,"28","","","","30","","","","","32"/)

		plot(0) 			= gsn_csm_contour_map(wks,dim_avg_n(temp(:,0,:,:),0),resT)


; ----------- PLOT SALINITY -------------------------------------
		;wks = gsn_open_wks ("x11","S")
		
		resS 				= res

  		;resS@tiMainString   		= "Salinity (10-3 PSU)" ; title
  		resS@lbTitleString 		= "Salinity (10~S~-3~N~ PSU)" ; bar title
		resS@gsnLeftString  		= " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ z="+dep(0)+"m ~C~ range=["+min(dim_avg_n(sws(:,0,:,:),0))+","+max(dim_avg_n(sws(:,0,:,:),0))+"]"; right title
  		resS@gsnRightString  		= "10~S~-3~N~ PSU"; right title
		
  		
		cmap 				= read_colormap_file("BkBlAqGrYeOrReViWh200")
		resS@cnFillPalette 		= cmap(18:177,:) ;ncolor=200
		
		;Select our own contour levels.
  		resS@cnLevelSelectionMode 	= "ManualLevels"
  		resS@cnLevelSpacingF      	= 0.4
		resS@cnMinLevelValF 		= 13.0
		resS@cnMaxLevelValF 		= 40.0

		plot(1) 			= gsn_csm_contour_map(wks,dim_avg_n(sws(:,0,:,:),0),resS)
		


; ----------- PLOT VELOCITY -------------------------------------

		;wks = gsn_open_wks ("x11","V")

		resV 				= res

		cmapV 				= read_colormap_file("WhiteBlueGreenYellowRed")
		resV@cnFillPalette 		= "WhiteBlueGreenYellowRed"

		; Title
		;resV@tiMainString 		= "Velocity (m/s) - feb-apr 1998"
  		resV@lbTitleString 		= "Current (m/s)" ; bar title
		resV@gsnLeftString  		= " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ z="+dep(0)+"m ~C~ range=["+min(dim_avg_n(spd(:,0,:,:),0))+","+max(dim_avg_n(spd(:,0,:,:),0))+"]"; right title
  		resV@gsnRightString  		= "m/s"; right title
		

		resV@mpGeophysicalLineColor     = "Grey18"
  		resV@mpGeophysicalLineThicknessF= 1.5	
	
		resV@vfXArray 			= lon
		resV@vfYArray 			= lat

		resV@cnLevelSelectionMode 	= "ManualLevels" ; set manual contour levels
		resV@cnMinLevelValF 		= 0.00
		resV@cnMaxLevelValF 		= 0.50
		resV@cnLevelSpacingF 		= 0.02
		;resV@lbLabelStrings 		= (/"0.00","","","","","0.20","","","",""\
	;					,"0.40","","","","","0.60","","","",""\
	;					,"0.80","","","","","1.00"/) 

		; vector specifications
		;resV@vcGlyphStyle 		= "CurlyVector"
		resV@vcRefMagnitudeF 		= 0.5
		resV@vcRefLengthF 		= 0.06
		resV@vcRefAnnoOrthogonalPosF 	= -1.0	; move ref vector
		resV@vcRefAnnoArrowLineColor 	= "black"
		resV@vcRefAnnoArrowUseVecColor 	= False	; don't use vec color for ref

		resV@vcLineArrowThicknessF 	= 0.1		; change vector thickness
		resV@vcVectorDrawOrder 		= "PostDraw"	; draw vectors at last
		resV@vcGlyphOpacityF 		= 0.7

		mean_u 				= dim_avg_n(u(:,0,:,:),0)
		mean_v 				= dim_avg_n(v(:,0,:,:),0)
		
		plot(2) 			= gsn_csm_vector_scalar_map_ce(wks,mean_u,mean_v\
							,dim_avg_n(spd(:,0,:,:),0),resV)
		

; --------------------- Draw panel ---------------------------------------
		
		gsn_panel(wks,plot,(/1,3/),False)




	end if

end
