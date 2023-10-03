
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% THIS CODE SIMULATES A SET OF simpleETAS SYNTHETIC EARTHQUAKE CATALOGS %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; 

%% SET PATH FOR SUPPORTING FUNCTIONS

addpath('/simplETAS-main/supporting_functions') ; 

%% SET DIRECTORIES

% absolute path of code working directory
wd = '/simplETAS-main/simulations' ;
% choose outputs directory
out_dir = '/simplETAS-main/simulations/outputs' ;

% absolute path of learning catalog
cat_dir = '/simplETAS-main/simulations/input_examples/HORUS_1972-2021_z30_M3+_no-etna_MPS19.txt' ;
% absolute path of mu(x,y) matrix (from national seismic hazard model)
bg_dir = '/simplETAS-main/simulations/input_examples/rate_norm_matrix_uncut.csv' ;
% absolute path of textfile with lat-lon vertices of authoritative region for mu(x,y)
% (in degrees, clockwise order starting from the top-left point) 
authreg_dir = '/simplETAS-main/simulations/input_examples/Polygon_Italy_rectangle_clockwise.dat' ;
% absolute path of time evolution
dt_dir = '/simplETAS-main/simulations/input_examples/dt_cumulative_392yrs.txt' ;

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

% load matrix representing the PDF of background events
bg_mtrx = readmatrix(bg_dir) ;

% load coordinates of authoritative region for mu(x,y)
auth_region = readmatrix(authreg_dir) ;

% load learning catalog of observations and trim it according to mag0:

% For "one shot" simulations, this catalog includes learning events occurred at t<=0.
% For retrospective and pseudo-prospective applications of the model to
% seismic sequences, this catalog includes pre-sequence (t<=0) + sequence-specific (t>0)seismicity
obs_raw = readmatrix(cat_dir) ;
obs = (obs_raw(obs_raw(:,4)>=mag0,:)) ; clear obs_raw ;
obsdim = size(obs) ;

% load vector indicating timing of model update since t=0 (in days)
dt = readmatrix(dt_dir) ;
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
        events_bg(j,9) = ceil(events_bg(j,1)) ; 
        events_bg(j,4) = round(events_bg(j,4),2) ;
    end
    
    parents_cat = cat(1,obs,events_bg) ;

    % then, simulate aftershocks
    for k = 1:b
        tmin_aft = dt(k) ;
        tmax_aft = dt(k+1) ;

        %parents_cat_filter = parents_cat(parents_cat(:,1)<=tmin_aft,1:8) ; % use tmin when simulations are updated with timesteps indicated in dt
        parents_cat_filter = parents_cat(parents_cat(:,1)<=tmax_aft,1:8) ;  % use tmax for "one shot" simulations
        
        parents_init = parents_cat_filter ;        % initial list of parent events
        parents_init = rmmissing(parents_init,1) ; % if any, remove NAs (then, check why you get them in the catalog ...)

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
    simulated_catalog = cat(1,events_aft, events_bg(:,1:8)) ;
    simulated_catalog = rmmissing(simulated_catalog,1) ;
    simulated_catalog = sortrows(simulated_catalog,1) ;

    % SAVE THE SIMULATION
    if i < 10
       writematrix(simulated_catalog,strcat(out_dir,'/Simulation_0000',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (10 <= i) && (i < 100) 
       writematrix(simulated_catalog,strcat(out_dir,'/Simulation_000',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (100 <= i) && (i < 1000)
       writematrix(simulated_catalog,strcat(out_dir,'/Simulation_00',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (1000 <= i) && (i <= 10000)
       writematrix(simulated_catalog,strcat(out_dir,'/Simulation_0',num2str(i),'.txt'),'Delimiter','tab') ; 
       elseif (10000 <= i) && (i <= 100000)
       writematrix(simulated_catalog,strcat(out_dir,'/Simulation_',num2str(i),'.txt'),'Delimiter','tab') ; 
    end

    % clear some variables
    clear j k events_bg bgdim parents_cat tmin_aft tmax_aft parents_cat_filter parents_init sim_filtered events_aft simulated_catalog

end

% Save the set of parameters in the same folder of the simulated catalogs
writematrix(theta,strcat(out_dir,'/Parameters.txt'),'Delimiter','tab') ; 

