load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	;data = addfile("data/dataset-global-reanalysis-phy-001-025-1997-98_dec-mar.nc","r")
	;data = addfile("data/dataset-global-reanalysis-phys-001-011-1997-98_dec-mar.nc","r")
	data = addfile("data/dataset-global-reanalysis-phys-001-011-1995_feb-apr.nc","r")
	
	; daily mean 2015-01-01 to 2015-02-01
	; temperature, salinity, velocity
	
	;print(data)

	i=1
	if (i.eq.1) then

; ----------- DATA -------------------------------------------------

		;sws = short2flt(data->salinity(:,0,:,:))	; sea water salinity
		;temp = short2flt(data->temperature(:,0,:,:)) 	; temperature
		;v = short2flt(data->v(:,0,:,:))		; northward velocity
		;u = short2flt(data->u(:,0,:,:))		; eastward velocity
		;lat = data->latitude				; latitude
		;lon = data->longitude				; longitude

		sws = data->sosaline		; sea water salinity
		temp = data->sosstsst 		; temperature
		v = data->somecrty		; northward velocity
		u = data->sozocrtx		; eastward velocity
		lat = data->lat			; latitude
		lon = data->lon			; longitude



		;temp = temp - 273.15	; K->C
		temp@units = "(C)"
                ;print(temp(:,:))
	

; ----------- PLOT TEMPERATURE -------------------------------------

		resources = True

		resources@cnMonoLineColor = False
		resources@tiMainString = "Temperature Angola region - feb-apr 1995"
		resources@cnFillOn = True
		resources@tiXAxisString = lon@long_name
		resources@tiYAxisString = lat@long_name
		resources@sfXArray = lon
		resources@sfYArray = lat
		resources@pmLabelBarDisplayMode = "Always"
		resources@lbPerimOn             = False
		resources@cnInfoLabelOn = False
	
		wks = gsn_open_wks ("png","angola_T_1995_feb-apr")

		plot = gsn_contour(wks,dim_avg_n(temp,0),resources)


; ----------- PLOT VELOCITY -------------------------------------
	
		resources = True
	
		resources@tiMainString = "Velocity Angola region - feb-apr 1995"
		resources@vcMinFracLengthF = 0.33
		resources@vcRefMagnitudeF  = 20.0
		resources@vcRefLengthF     = 0.045
		resources@vcMonoLineArrowColor  = False   ; Draw vectors in color.
		resources@tiXAxisString = lon@long_name
		resources@tiYAxisString = lat@long_name

		resources@vfXCStartV  = lon(0,0)               ; Define X/Y axes range
		resources@vfXCEndV    = lon(dimsizes(lon)-1,0) ; for vector plot.
		resources@vfYCStartV  = lat(0,0)
		resources@vfYCEndV    = lat(dimsizes(lat)-1,0)

		resources@pmLabelBarDisplayMode    = "Always"
		resources@lbPerimOn                = False
		resources@vcGlyphStyle = "CurlyVector"
		resources@cnInfoLabelOn = False

		wks = gsn_open_wks ("png","angola_V_1995_feb-apr")

		plot = gsn_vector(wks,dim_avg_n(u,0),dim_avg_n(v,0),resources)



		delete(plot)	; Clean up
		delete(temp)
		delete(resources)

	end if

end
