load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

	data = addfile("data/dataset-global-reanalysis-phy-001-025-01031995_deep.nc","r")
	
	; temperature, salinity, velocity on a section plot
	
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

		lat_section = -2
		lat_ind = ind(lat.eq.lat_section)
		if (dimsizes(dep).lt.2)
			print("***Data must have depth***")
			exit
		end if
		

; ----------- INITIALIZATION General resourses -----------------------------------

		res 				= True

		plot 				= new(4,graphic)
		file_title 			= ""+dd(0)+"-"+mm(0)+"-"+yyyy(0)+"_"+dd(dimsizes(dd)-1)+"-"+mm(dimsizes(dd)-1)+"-"+yyyy(dimsizes(dd)-1)+"_lat_"+lat_section+""
		wks = gsn_open_wks("ps",file_title)

		res@gsnDraw  			= False			; don't draw
		res@gsnFrame 			= False			; don't advance frame
		res@gsnAddCyclic          	= False       		; data already has cyclic point
		;res@gsnMaximize 		= "True"		; Maximize plot in frame
		res@gsnLeftStringFontHeightF 	= 0.015
		res@gsnLeftStringOrthogonalPosF = 0.002
		res@gsnRightStringFontHeightF 	= 0.015
		res@gsnRightStringOrthogonalPosF= 0.006
  		;res@gsnScalarContour	 	= True			; contours desired

		res@cnMonoLineColor 		= False
		res@cnFillOn 			= True			; Fill with color
		res@cnLinesOn 			= False			; No contour lines
		res@cnInfoLabelOn 		= False
		
		res@pmTickMarkDisplayMode 	= "Always" 		; Default lat/lon labels
	
		res@tiMainOffsetYF 		= -0.03
		res@tiYAxisString 		= "Depth (m)"
		res@tiXAxisString 		= "Longitude"	

		res@lbOrientation        	= "Vertical"
		res@lbTitlePosition      	= "Right"               ; title location
  		res@lbTitleDirection     	= "Across"              ; letter angle
  		res@lbTitleAngleF        	= 90.                   ; title angle
  		res@lbTitleFontHeightF   	= 0.02                  ; font height
		res@lbTitleOffsetF	 	= 0.3
		res@lbPerimOn            	= False
		
		res@sfXArray 			= lon
		res@sfYArray 			= dep
		res@trYReverse 			= True             	; reverses y-axis

		;res@pmLabelBarDisplayMode 	= "Always"
		res@tmXBLabelFontHeightF 	= 0.018
		res@tmYLLabelFontHeightF 	= 0.018

; ----------- PLOT TEMPERATURE -------------------------------------

		;wks = gsn_open_wks ("x11","T")

		resT = res		
		
		cmap = read_colormap_file("BkBlAqGrYeOrReViWh200")
		resT@cnFillPalette = cmap(18:177,:) ;ncolor=200

		;resT@tiMainString = "Sea Surface Temperature (C)"
		resT@gsnLeftString  = " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ lat="+lat_section+"~S~o~N~ ~C~ range=["+min(dim_avg_n(temp(:,:,lat_ind,:),0))+","+max(dim_avg_n(temp(:,:,lat_ind,:),0))+"]"; right title
  		resT@gsnRightString  = "~S~o~N~C"; right title
		resT@lbTitleString = "Temperature (~S~o~N~C)" ; bar title
		
		;Select our own contour levels.
  		resT@cnLevelSelectionMode = "ManualLevels"
  		resT@cnLevelSpacingF      = 0.4
		resT@cnMinLevelValF = 12.0
		resT@cnMaxLevelValF = 32.4

		plot(0) = gsn_csm_contour_map(wks,dim_avg_n(temp(:,:,lat_ind,:),0),resT)



; ----------- PLOT SALINITY -------------------------------------
		;wks = gsn_open_wks ("x11","S")
		
		resS = res

  		;resS@tiMainString   = "Salinity (10-3 PSU)" ; title
  		resS@lbTitleString = "Salinity (10~S~-3~N~ PSU)" ; bar title
		resS@gsnLeftString  = " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ lat="+lat_section+"~S~o~N~ ~C~ range=["+min(dim_avg_n(sws(:,:,lat_ind,:),0))+","+max(dim_avg_n(sws(:,:,lat_ind,:),0))+"]"; right title
  		resS@gsnRightString  = "10~S~-3~N~ PSU"; right title
		
		cmap = read_colormap_file("BkBlAqGrYeOrReViWh200")
		resS@cnFillPalette = cmap(18:177,:) ;ncolor=200
		;resS@vcLevelPalette = cmap(18:177,:)

		;Select our own contour levels.
  		resS@cnLevelSelectionMode = "ManualLevels"
  		resS@cnLevelSpacingF      = 0.4
		resS@cnMinLevelValF = 13.0
		resS@cnMaxLevelValF = 40.0 

		plot(1) = gsn_csm_contour(wks,dim_avg_n(sws(:,:,lat_ind,:),0),resS)


; ----------- PLOT VELOCITY -------------------------------------

		;wks = gsn_open_wks ("x11","V")

		resV = res

		cmapV = read_colormap_file("BlRe")
		resV@cnFillPalette = "BlRe"

		; Title
		;resV@tiMainString = "Velocity (m/s) - feb-apr 1998"
  		resV@gsnRightString  = "m/s"; right title

		resV@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
		resV@cnMinLevelValF = -0.50
		resV@cnMaxLevelValF = 0.50
		resV@cnLevelSpacingF = 0.02

		mean_u = dim_avg_n(u(:,:,lat_ind,:),0)
		mean_v = dim_avg_n(v(:,:,lat_ind,:),0)
		

		resV@lbTitleString = "Zonal velocity (m/s)" ; bar title
		resV@gsnLeftString  = " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ lat="+lat_section+"~S~o~N~ ~C~ range=["+min(dim_avg_n(u(:,:,lat_ind,:),0))+","+max(dim_avg_n(u(:,:,lat_ind,:),0))+"]"; right title
		plot(2) = gsn_csm_contour(wks,mean_u,resV)

		delete(resV@lbTitleString)
		delete(resV@gsnLeftString)

		resV@lbTitleString = "Meridional velocity (m/s)" ; bar title
		resV@gsnLeftString  = " Mean from "+mm(0)+"/"+dd(0)+"/"+yyyy(0)+" to "+mm(dimsizes(mm)-1)+"/"+dd(dimsizes(mm)-1)+"/"+yyyy(dimsizes(mm)-1)+" ~C~ lat="+lat_section+"~S~o~N~ ~C~ range=["+min(dim_avg_n(v(:,:,lat_ind,:),0))+","+max(dim_avg_n(v(:,:,lat_ind,:),0))+"]"; right title
		plot(3) = gsn_csm_contour(wks,mean_v,resV)
		

; --------------------- Draw panel ---------------------------------------
		
		gsn_panel(wks,plot,(/4,1/),False)




	end if

end
