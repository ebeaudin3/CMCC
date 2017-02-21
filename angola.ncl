load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	;global_data = addfile("data/dataset-global-reanalysis-phy-001-025-ran-fr-glorys2v4-daily_1486999284990.nc","r")
	global_data = addfile("data/dataset-global-reanalysis-phy-001-025-ran-fr-glorys2v4-daily_1487083789645.nc","r")
	; daily mean 2015-01-01 to 2015-02-01
	; temperature, salinity, velocity
	
	;print(global_data)

	i=1
	if (i.eq.1) then

; ----------- DATA -------------------------------------------------

		dep_ind = 0

		mlp = short2flt(global_data->mlp)		; mixed layer thickness
		sws = short2flt(global_data->salinity)		; sea water salinity
		ssh = short2flt(global_data->ssh)		; sea surface height
		dep = global_data->depth			; depth
		temp = short2flt(global_data->temperature(0,dep_ind,:,:)) ; temperature
		v = short2flt(global_data->v(0,dep_ind,:,:))		; northward velocity
		u = short2flt(global_data->u(0,dep_ind,:,:))		; eastward velocity
		time = global_data->time			; time
		lat = global_data->u&latitude			; latitude
		lon = global_data->v&longitude			; longitude

		temp = temp - 273.15	; K->C
		temp@units = "(C)"
                ;print(temp(:,:))
	

; ----------- PLOT TEMPERATURE -------------------------------------

		resources = True

		resources@cnMonoLineColor = False
		resources@tiMainString = "Temperature Angola region - daily average Jan 1st, 2015"
		resources@cnFillOn = True
		resources@tiXAxisString = lon@long_name
		resources@tiYAxisString = lat@long_name
		resources@sfXArray = lon
		resources@sfYArray = lat
		resources@pmLabelBarDisplayMode = "Always"
		resources@lbPerimOn             = False
		resources@cnInfoLabelOn = False
	
		wks = gsn_open_wks ("x11","angola_temp_big")

		plot = gsn_contour(wks,temp,resources)


; ----------- PLOT VELOCITY -------------------------------------
	
		resources = True
	
		resources@tiMainString = "Velocity Angola region - daily average Jan 1st, 2015"
		resources@vcMinFracLengthF = 0.33
		resources@vcRefMagnitudeF  = 20.0
		resources@vcRefLengthF     = 0.045
		resources@vcMonoLineArrowColor  = False   ; Draw vectors in color.
		resources@tiXAxisString = lon@long_name
		resources@tiYAxisString = lat@long_name

		resources@vfXCStartV  = lon(0)               ; Define X/Y axes range
		resources@vfXCEndV    = lon(dimsizes(lon)-1) ; for vector plot.
		resources@vfYCStartV  = lat(0)
		resources@vfYCEndV    = lat(dimsizes(lat)-1)

		resources@pmLabelBarDisplayMode    = "Always"
		resources@lbPerimOn                = False
		resources@vcGlyphStyle = "CurlyVector"
		resources@cnInfoLabelOn = False

		wks = gsn_open_wks ("x11","angola_vel_big")

		plot = gsn_vector(wks,u,v,resources)



		delete(plot)	; Clean up
		delete(temp)
		delete(resources)

	end if

end
