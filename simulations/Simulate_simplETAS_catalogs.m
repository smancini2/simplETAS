%% From Mancini and Marzocchi (2023), SRL

%% Version of Dec 12, 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% THIS CODE SIMULATES A SET OF simplETAS SYNTHETIC EARTHQUAKE CATALOGS %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE CODE READS EARTHQUAKES CATALOGS IN ZMAP-like FORMAT

% i.e., the learning earthquake catalog must have the following format (columns):
% column 1: longitude [deg]
% column 2: latitude [deg]
% column 3: year
% column 4: month
% column 5: day
% column 6: magnitude
% column 7: depth [km]
% column 8: hour
% column 9: minute
% column 10: second

% example: [13.2335 42.6983 2016 08 24 6.18 8.10 01 36 32]

% THE OUTPUT SIMULATED CATALOGS WILL BE IN "ENRICHED" ZMAP FORMAT

% i.e., featuring the following additional columns
% column 11: unique event index 
% column 12: type of simulated event (1=triggered, 2=background)
% column 13: event generation (if =0, indicates a background event)
% column 14: index of the related parent event (if =0, indicates a background event)

% NOTE: these are 2D simulations, so the depth (column 7) will be set=0km


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; 

%% SET PATH TO THE CODE FOLDER AND SUBFOLDERS

addpath(genpath('/simplETAS-main')) ; % set your full path to 'simplETAS-main' folder 

%% SET DIRECTORIES

% absolute path of code working directory
wd = '/simplETAS-main/simulations' ;
% choose outputs directory
out_dir = '/simplETAS-main/simulations/outputs' ;

% absolute path of learning catalog
cat_dir = '/simplETAS-main/simulations/input_examples/HORUS_1972-2021_z30_M3+_no-etna_MPS19_ZMAP.txt' ;
% absolute path of mu(x,y) matrix (from national seismic hazard model)
bg_dir = '/simplETAS-main/simulations/input_examples/rate_norm_matrix_uncut.csv' ;
% absolute path of textfile with lat-lon vertices of authoritative region for mu(x,y)
% (in degrees, clockwise order starting from the top-left point) 
authreg_dir = '/simplETAS-main/simulations/input_examples/Polygon_Italy_rectangle_clockwise.dat' ;

cd(wd) ;

%% SET MODEL INPUT PARAMETERS

% GENERIC PARAMETERS FOR SIMULATIONS
nsim = 1000 ;        % number of simulations
K = 0 ;              % initialize the variable "K" used in Simulate_aftershocks.m

reflat = 38.7328 ;   % used for to express x-y [in km] in learning catalog and in all simulation outputs
reflon = 13.1237 ;   % used for to express x-y [in km] in learning catalog and in all simulation outputs
binning = 0.1 ;      % spatial binning in degrees, square bins required
xbin = lldistkm([reflat reflon],[reflat reflon+binning]) ; % get actual spatial binning in km along longitude
ybin = lldistkm([reflat reflon],[reflat+binning reflon]) ; % get actual spatial binning in km along latitude

fore_start = [2022 01 01 00 00 00] ;  % start of the forecast window [yyyy mm dd hh min sec] 
fore_duration = 143178 ;              % [in days] total forecast horizon 
fore_upd_freq = 143178 ;              % [in days] forecast update frequency -- for "one shot" forecasts, fore_duration == fore_steps

% MAGNITUDE PARAMETERS
bvalue = 1.0 ;  % assumed b-value in GR distribution
mbin = 0.1 ;    % magnitude bin size
mag0 = 3.95 ;   % minimum triggering magnitude

% limits of the truncated GR distribution for the magnitudes of the simulated events 
% (i.e., the magnitude range of events to be simulated)
mag_large = 7.5 ;  % maximum magnitude to simulate
mag_small = 4.0 ;  % input value for magnitude simulation --> the minumum magnitude to simulate is = mag_small-(mbin/2)
                   % for instance, if the minimum magnitude to be simulated is 3.95, then mag_small = 4.0 and mbin = 0.1

% simplETAS PARAMETERS
v = 18.2713 ; % in events/year
A = 0.0473 ; 
alpha = bvalue*log(10) ; beta = bvalue*log(10) ; 
p = 1.15 ; 
c = 0.005 ; % in days 
gamma = 1.5 ; q = 1.5 ; 
D = 1 ; % in km^2

% create matrices with the set of simplETAS parameters for the simulations
omega = [beta,mag_small] ;
theta = [A,c,alpha,p,gamma,D,q,beta,mag_small] ;


%% OTHER PRELIMINARY STEPS

% Load learning catalog of observations and change format (for internal use only)

% NOTE: For "one shot" simulations, this catalog includes learning events occurred at t<=fore_start only.
% For retrospective and pseudo-prospective applications of the model to seismic sequences, 
% this catalog includes pre-sequence (t<=fore_start) + sequence-specific (t>fore_start) seismicity

obs_raw_ZMAP = readmatrix(cat_dir) ;
obsdim = size(obs_raw_ZMAP) ;

for n = 1:obsdim(1)
    obs_raw(n,1) = [etime([obs_raw_ZMAP(n,[3 4 5 8 9 10])],fore_start)]/(24*60*60) ; %y mon day hr min sec
end

[obs_raw(:,3),obs_raw(:,2)] = geo2merc(obs_raw_ZMAP(:,2),obs_raw_ZMAP(:,1),reflon,reflat) ;
obs_raw(:,4) = obs_raw_ZMAP(:,6) ; 
obs_raw(:,5) = 2 ; % this means that all entries are background events
obs_raw(:,6) = 0 ; % put generation = 0 for all entries
for n = 1:obsdim(1)
    obs_raw(n,7) = n ; 
end
obs_raw(:,8) = 0 ; 

clear obsdim n obs_raw_ZMAP

% load matrix representing the PDF of background events
bg_mtrx = readmatrix(bg_dir) ;

% load coordinates of authoritative region for mu(x,y)
auth_region = readmatrix(authreg_dir) ;

% Cut the learning catalog according to mag0:
obs = (obs_raw(obs_raw(:,4)>=mag0,:)) ; clear obs_raw ;

% create vector indicating the regular timing of model update since fore_start [in days, incuding t=0]

dt = (0:fore_upd_freq:fore_duration)' ;dtdim = size(dt) ;
dtdim = size(dt) ;
b = dtdim(1)-1 ; % this line corresponds to the desired total number of forecast windows


%% SIMULATE simplETAS CATALOGS

avg_bg = (dt(end,1)/365.25)*v ; % average number of background events to simulate (lambda parameter for bg simulation)
tmin_bg = dt(1,1) ;
tmax_bg = dt(end,1) ;

for i = 1:nsim

    % first, simulate a catalog of background events
    events_bg = Simulate_background(avg_bg, bg_mtrx, auth_region, xbin, ybin, mbin, ...
        reflat, reflon, omega, mag_large, tmin_bg, tmax_bg) ;

    % bind together the catalog of learning events (t<=0) and that of background
    % events just simulated (t>0), to create a final catalog of "parent"
    % events
    events_bg = events_bg(events_bg(:,4)>=mag0,:) ; 
    bgdim = size(events_bg) ;
    
    for j = 1:bgdim(1)
        events_bg(j,7) = max(obs(:,7))+j ;
        %events_bg(j,9) = ceil(events_bg(j,1)) ;  
        events_bg(j,4) = round(events_bg(j,4),2) ;
    end
    
    parents_cat = cat(1,obs,events_bg) ;

    % then, simulate aftershocks
    for k = 1:b
        tmin_aft = dt(k) ;
        tmax_aft = dt(k+1) ;

        if dtdim(1) == 2   % (for "one shot" simulations, i.e. with one tmin and one tmax only)
           parents_cat_filter = parents_cat(parents_cat(:,1)<=tmax_aft,1:8) ;
        
        else               % (i.e., for simulations that are updated a number of times using the timesteps indicated in dt)
           parents_cat_filter = parents_cat(parents_cat(:,1)<=tmin_aft,1:8) ;
        end
        
        parents_init = parents_cat_filter ;         % initial list of parent events
        parents_init = rmmissing(parents_init,1) ;  % if any, remove NAs (then, check why you get them in the catalog ...)

        fin_gen = Simulate_aftershocks(parents_init,theta,mag_large,mbin,tmin_aft,tmax_aft) ; % do simulation
        
        sim_filtered = fin_gen(fin_gen(:,8)>0,:) ;
        clear fin_gen ;

        if k == 1
           events_aft = sim_filtered ;
        else
           events_aft = cat(1,events_aft, sim_filtered) ;
        end
 
    end
    
    % bind simulated background events + simulated aftershocks of all
    % generations, to get the final (total) simulated earthquake catalog
    simulated_catalog = cat(1,events_aft, events_bg) ;
    simulated_catalog = rmmissing(simulated_catalog,1) ;
    simulated_catalog = sortrows(simulated_catalog,1) ;

    % back to "enriched" ZMAP format 
    [lat,lon] = merc2geo(simulated_catalog(:,3),simulated_catalog(:,2),reflon,reflat) ;
    lat(lat>90)=90 ; lat(lat<-90)=-90 ;
    lon(lon>180)=180 ; lon(lon<-180)=-180 ;
    temp(:,2) = lat ; temp(:,1) = lon ;
    [temp(:,3) temp(:,4) temp(:,5) temp(:,8) temp(:,9) temp(:,10)] = datevec(datenum(fore_start)+simulated_catalog(:,1)) ; 
    temp(:,6) = simulated_catalog(:,4) ;
    temp(:,7) = 0 ;
    temp(:,11) = simulated_catalog(:,7) ;
    temp(:,12) = simulated_catalog(:,5) ;
    temp(:,13) = simulated_catalog(:,6) ;
    temp(:,14) = simulated_catalog(:,8) ;

    simulation_final = temp ;
    clear temp simulated_catalog

    % SAVE THE SIMULATION
    if i < 10
       writematrix(simulation_final,strcat(out_dir,'/Simulation_0000',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (10 <= i) && (i < 100) 
       writematrix(simulation_final,strcat(out_dir,'/Simulation_000',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (100 <= i) && (i < 1000)
       writematrix(simulation_final,strcat(out_dir,'/Simulation_00',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (1000 <= i) && (i <= 10000)
       writematrix(simulation_final,strcat(out_dir,'/Simulation_0',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (10000 <= i) && (i <= 100000)
       writematrix(simulation_final,strcat(out_dir,'/Simulation_',num2str(i),'.txt'),'Delimiter','tab') ; 
    end

    % clear some variables
    clear j k events_bg bgdim parents_cat tmin_aft tmax_aft parents_cat_filter parents_init sim_filtered events_aft simulation_final n lat lon

end

% Save the set of parameters in the same folder of the simulated catalogs
writematrix(theta,strcat(out_dir,'/Parameters.txt'),'Delimiter','tab') ; 

