function [ylat,xlon]=geo2merc(lat,lon,xcen,ycen)

%Finding the distance between 2 lat/lon points
ylat=(lat-ycen)*110.94;

xlon=(lon-xcen)*(110.94*cosd(ycen));

















