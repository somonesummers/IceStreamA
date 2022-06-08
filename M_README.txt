11 column CSV
M = [Xi(:),Yi(:),lat(:),lon(:),u(:),v(:),spd(:),b_raw(:),sf_raw(:),thickness(:),ice_mask(:)];
(0) Xi = NSIDC Polar Stereographic Grid northing (https://nsidc.org/data/polar-stereo/ps_grids.html) [meters]
(1) Yi = NSIDC Polar Stereographic Grid easting (https://nsidc.org/data/polar-stereo/ps_grids.html) [meters]
(2) Latitude
(3) Longitude
(4) ice surface x velocity (easting) [m/yr]
(5) ice surface y velocity (northing) [m/yr]
(6) ice surface speed m[/yr]
(7) bed elevation (m) relative to gl04c geoid [m]
(8) ice surface elevation (m) relative to gl04c geoid [m]
(9) ice thickness [m]
(10) ice mask 0 = grounded, 1 = ice shelf, 127 = ocean 