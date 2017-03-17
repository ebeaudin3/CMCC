load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

	ncl_data = addfile("data/bathymetry/GEBCO_2014_2D_2.0_-19.0_16.0_2.0.nc","r")

	;print(ncl_data)

	elev = ncl_data->elevation	; elevation
	lat = ncl_data->lat		; latitude
	lon = ncl_data->lon		; longitude

; ----------- INITIALIZATION -------------------------------------

	res = True

	res@gsnPanelLabelBar = True
	res@lbOrientation = "Vertical"
	res@cnMonoLineColor = False
	res@cnFillOn = True
	res@cnLinesOn = False

	res@pmTickMarkDisplayMode = "Always" 	; Default lat/lon labels
	;res@gsnMaximize = "True"		; Maximize plot in frame
	res@gsnAddCyclic          = False       ; data already has cyclic point

	res@tiMainOffsetYF = -0.03
	res@gsnRightStringOrthogonalPosF = -0.002
  	res@gsnLeftString  = ""              ; left title
	res@gsnRightString  = "Bathymetry" ; right title

	; Zoom in on map area that we are interested in
	res@mpLimitMode      = "Corners"
  	res@mpLeftCornerLatF = lat(0)
  	res@mpLeftCornerLonF = lon(0)
  	res@mpRightCornerLatF = lat(dimsizes(lat)-1)
  	res@mpRightCornerLonF = lon(dimsizes(lon)-1)
	res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
	res@mpDataBaseVersion = "HighRes"

	res@lbOrientation        = "Vertical"
	res@lbTitleString 	 = "Depth (m)" ; bar title
	res@lbTitlePosition      = "Right"                           ; title location
  	res@lbTitleDirection     = "Across"                          ; letter angle
  	res@lbTitleAngleF        = 90.                               ; title angle
  	res@lbTitleFontHeightF   = 0.02                              ; font height
	res@lbTitleOffsetF	 = 0.3


; ----------- PLOT BATHYMETRY -------------------------------------
	
	wks = gsn_open_wks ("x11","angola_bathymetry")

	ccmap = read_colormap_file("GMT_drywet")
	cmap = ccmap(::-1,:)
	res@cnFillPalette = cmap
	res@vcLevelPalette = cmap

	res@cnMonoLineColor = False
	res@tiMainString = "Coastal Angola Region"
	;res@tiXAxisString = lon@long_name
	;res@tiYAxisString = lat@long_name
	res@sfXArray = lon
	res@sfYArray = lat
	res@pmLabelBarDisplayMode = "Always"
	res@lbPerimOn             = False
	
	;Select our own contour levels.

  	res@cnLevelSelectionMode = "ManualLevels"
  	res@cnLevelSpacingF      = 200

	base = gsn_csm_contour_map(wks,elev,res) 
   	draw(base)  

	;plot = gsn_contour(wks,elev,res)

end
