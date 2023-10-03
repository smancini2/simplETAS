%% From Mancini and Marzocchi (2023), SRL

%% Version of Aug 28, 2023
%
% Acknowledgment: some functions and the main structure of this code were originally written as 
% part of an R code by Dr Stefanie Seif and used for the following paper:

% Seif, S., A. Mignan, J. D. Zechar, M. J. Werner, and S. Wiemer (2017). 
% Estimating ETAS: The effects of truncation, missing data, and model assumptions, 
% J. Geophys. Res. Solid Earth 121, 449-469. https://doi.org/10.1002/2016JB012809.  

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
% THIS FUNCTION SIMULATES simplETAS BACKGROUND EVENTS %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function events_bg = simulate_background(avg_bg, bg_mtrx, auth_region, xbin, ybin, mbin, ...
         reflat, reflon, omega, mag_large, tmin_bg, tmax_bg)

%bg_mtrx = flip(bg_mtrx,1) ;
bg_mtrx = transpose(flip(bg_mtrx,1)) ;

% Target Region [anti clockwise]
px = flip(auth_region(:,1)) ;
py = flip(auth_region(:,2)) ;

% Append first element as the last --> necessary for integration!
px = cat(1,px,px(1)) ;
py = cat(1,py,py(1)) ;

% Express px and py relative to the center of the region
px = (px-reflon)*Lon2Km(reflat) ;
py = (py-reflat)*6378.1*pi/180 ; 

% Get grid extensions
[xpol(1),xpol(2)] = bounds(px) ; % in km relative to the center
[ypol(1),ypol(2)] = bounds(py) ; 

xpol_round = xpol ;
ypol_round = ypol ;
%x.pol.round = c(floor(x.pol[1]/grid.bin)*grid.bin,ceiling(x.pol[2]/grid.bin)*grid.bin)
%y.pol.round = c(floor(y.pol[1]/grid.bin)*grid.bin,ceiling(y.pol[2]/grid.bin)*grid.bin)

% Divide the grid
xgrid = linspace(xpol_round(1),xpol_round(2),size(bg_mtrx,1))' ;
ygrid = linspace(ypol_round(1),ypol_round(2),size(bg_mtrx,2))' ;

xvec = repmat(xgrid,size(bg_mtrx,2),1) ;
yvec = repelem(ygrid,size(bg_mtrx,1)) ;

%bg_tar = nm2N(bg_mtrx) ;
bg_tar = nm2N(bg_mtrx') ; % only if line 12 is uncommented

% Transform rate to density
dens_mtrx = bg_tar/sum(bg_tar) ;
dens_cum = N2nm(cumsum(dens_mtrx),size(bg_mtrx,2),size(bg_mtrx,1)) ;


% Generate background
nbg = poissrnd(avg_bg) ;
tbg = tmin_bg+(tmax_bg-tmin_bg)*rand(nbg,1) ;
mbg = GRTruncSimul(omega(1),mag_large,omega(2),mbin, nbg) ;

for j = 1:nbg
    x = rand ;
    [row,col] = find(dens_cum>=x);
    GridInd(1,j) = row(1) ; GridInd(2,j) = col(1) ;
    clear row col x
end

xlist = xgrid(GridInd(1,:)) + xbin*rand(size(GridInd,2),1) ; %[km]
ylist = ygrid(GridInd(2,:)) + ybin*rand(size(GridInd,2),1) ; %[km]

temp = (1:1:nbg)' ;

events_bg(:,1:8) = [tbg,xlist,ylist,mbg,repmat(2,nbg,1),zeros(nbg,1),temp,zeros(nbg,1)] ;
events_bg = rmmissing(events_bg,1) ;

clear temp

end
