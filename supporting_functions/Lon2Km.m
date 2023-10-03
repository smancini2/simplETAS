
function out = Lon2Km(lat)

    rad_earth = 6378.1 ;  % [km]
    out = (rad_earth*cos(lat*pi/180)*pi/180) ;