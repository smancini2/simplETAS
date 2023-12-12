%% From Mancini and Marzocchi (2023), SRL

%% Version of Dec 12, 2023
%
% Other papers cited in the code
% E. Lippiello,  F. Giacco, L. de Arcangelis, W. Marzocchi, C. Godano (2014). 
% Parameter estimation in branching processes: Approximations and novel methods. 
% Bull. Seismol. Soc. Am., 104, 985-994.
%
% Schoenberg, F. P. (2013). Facilitated estimation of ETAS, 
% Bull. Seismol. Soc. Am. 103, no. 1, 601–605, doi: 10.1785/0120120146.

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% THIS CODE ESTIMATES THE simplETAS PARAMETERS AND THEIR UNCERTAINTY    %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
close all; clear all

%% SET PATH FOR SUPPORTING FUNCTIONS

addpath(genpath('/simplETAS-main')) ; % set your full path to 'simplETAS-main' folder 

%% 
 
% INPUT:  1) earthquake catalog in ZMAP format [obs_raw, 10 columns:]

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
     
%         2) file of the region coordinates [forecasting_region; 2 columns:
%                                          i)  longitude in degrees ; 
%                                          ii) latitude in degrees]
% 
% OUTPUT: 1) the earthquake catalog inside the 
%            space-time window [obs_raw_inside, 4 columns: 
%               i)   time in days since the first event in the whole catalog; 
%               ii)  longitude in degrees; 
%               iii) latitude in degrees; 
%               iv)  magnitude]
%         2) the estimated simplETAS free parameters

%% SET DIRECTORIES

% absolute path of the seismic catalog
cat_dir = '/simplETAS-main/parameters_estimation/input_examples/HORUS_1972-2021_z30_no-etna_ZMAP.txt' ; 

% absolute path of textfile with lon-lat vertices of authoritative region
% (in degrees, clockwise order starting from the top-left point)
forereg_dir = '/simplETAS-main/parameters_estimation/input_examples/Polygon_Italy_MPS19_clockwise.dat' ;

% absolute path of mu(x,y) matrix (from national seismic hazard model)
bg_dir = '/Users/simone/Desktop/simplETAS-main/parameters_estimation/input_examples/rate_norm.txt' ; 


%% LOAD INPUT CATALOG AND SELECT THE EARTHQUAKES IN THE SPACE-TIME WINDOW OF INTEREST

reflat = 38.7328 ;   % used for to transform x-y [in km] in learning catalog
reflon = 13.1237 ;   % used for to transform x-y [in km] in learning catalog

forecasting_region = readmatrix(forereg_dir) ;   

obs_raw_ZMAP = readmatrix(cat_dir) ;
obsdim = size(obs_raw_ZMAP) ;
last_event = obs_raw_ZMAP(end,[3 4 5 8 9 10]) ;

for n = 1:obsdim(1)
    obs_raw(n,1) = [etime([obs_raw_ZMAP(n,[3 4 5 8 9 10])],last_event)]/(24*60*60) ; %y mon day hr min sec
end
obs_raw(:,2) = obs_raw_ZMAP(:,1) ;
obs_raw(:,3) = obs_raw_ZMAP(:,2) ;
obs_raw(:,4) = obs_raw_ZMAP(:,6) ;
clear obs_raw_ZMAP n

%t0 = 0 ;
Nyears = 3 ;           % start of calculation for the likelihood (i.e., diration of learning period)
t1 = Nyears*365 ;      % the period for the estimation starts Nyears after the start of the catalog
t2 = -obs_raw(1,1) ;   % end of estimation period

inside_raw = inpolygon(obs_raw(:,2),obs_raw(:,3),forecasting_region(:,1),forecasting_region(:,2)) ; % Determining the earthquakes of the catalog inside the forecasting region
obs_raw_inside1 = obs_raw((inside_raw),:) ;    % inside the spatial region
obs_raw_inside = obs_raw_inside1(find((obs_raw_inside1(:,1)-obs_raw(1,1))>t1),:) ;      % inside the temporal window

%% SET MINIMUM MAGNITUDE mag0 AND CUT THE CATALOG ACCORDINGLY 

mag0 = 3.95 ;   % minimum triggering magnitude
obs = (obs_raw(obs_raw(:,4) >= mag0,:)) ; % cut the catalog according to mag0

time = obs(:,1)-obs_raw(1,1) ;
obs_inside = (obs_raw_inside(obs_raw_inside(:,4) >= mag0,:)) ; % cut the inside catalog according to mag0:
time_inside = obs_inside(:,1)-obs_raw(1,1) ;

N_eqk = length(obs) ;
N_eqk_inside = length(obs_inside) ;

clear obs_raw obs_raw_inside1 obs_raw_inside ;


%% SET THE BACKGROUND (from national seismic hazard model)
%  OUTPUT: bg_mtrx_inside, that is the background (bg_mtrx) inside the space-time window. 3 columns: i) longitude, ii) latitude; iii) normalized spatial probability

%  It may be very similar to the input background bg_mtrx if the latter is calculated and normalized on the same space window

bg_mtrx = readmatrix(bg_dir) ;
inside_bg = inpolygon(bg_mtrx(:,1),bg_mtrx(:,2),forecasting_region(:,1),forecasting_region(:,2)) ;
bg_mtrx_inside = bg_mtrx(inside_bg,:) ;


%% simplETAS FIXED PARAMETERS

bvalue = 1.0 ;  % assumed b-value in GR distribution
alpha = bvalue*log(10) ; beta = bvalue*log(10) ; 
p = 1.15 ; 
c = 0.005 ; % in days 
gamma = 1.5 ; q = 1.5 ; 
D = 1 ; % in km^2


%% CALCULATING SOME QUANTITIES USED FOR THE EQUATION TO BE SOLVED TO FIND THE PARAMETERS 
% (see equations in the supplementary material)

% Calculating H
%h = exp(alpha*(obs(:,4)-mag0)) ;    % approximation by Schoenberg (2013)
h = exp(alpha*(obs(:,4)-mag0)).*(1-((t2-time)/c+1).^(1-p)) ;   % this contains the approximation on the space integral as in Lippiello et al. BSSA 2014
%h = exp(alpha*(obs(:,4)-mag0)).*(1-((t2-time)/c+1).^-p) ;     % Error on equation 14 in Lippiello et al's (2014) paper
H = sum(h) ;

% Calculating Delta
du = bg_mtrx_inside(:,3) ;
Delta = sum(du) ;

% Calculating C(j)
for j = 1:N_eqk_inside
    cdum = 0 ;
    [ykm_j xkm_j] = geo2merc(obs_inside(j,3),obs_inside(j,2),reflon,reflat) ;
    for i = 1:N_eqk
        if time (i) < time_inside(j)  
            [ykm_i xkm_i] = geo2merc(obs(i,3),obs(i,2),reflon,reflat) ;
           omori = (c^(p-1)*(p-1))/(time_inside(j)-time(i)+c)^p ;
           utsu = exp(alpha*(obs(i,4)-mag0)) ;
           dist_decay = (q-1)/(pi*D*exp(gamma*(obs(i,4)-mag0)))* ...
               (1+((xkm_i-xkm_j)^2+(ykm_i-ykm_j)^2)/...
               (D*exp(gamma*(obs(i,4)-mag0))))^(-q) ;
           cdum = cdum+omori*utsu*dist_decay ; 
        end
    end
    C(j) = cdum ;
end

disp('Step 1')

% Calculating mu(j)

[bgkm_y bgkm_x] = geo2merc(bg_mtrx_inside(:,2),bg_mtrx_inside(:,1),reflon,reflat) ;
[eqkm_y eqkm_x] = geo2merc(obs_inside(:,3),obs_inside(:,2),reflon,reflat) ;
XX = [bgkm_x bgkm_y] ;
YY = [eqkm_x eqkm_y] ;
Idx = knnsearch(XX,YY) ;
mu = bg_mtrx_inside(Idx,3) ;

%
% Checked mu with this procedure. Results are the same.
%
% wgs84km = wgs84Ellipsoid("kilometer")
% for j=1:N_eqk_inside
%     if floor(j/100)*100 == j 
%         disp('iter=') 
%         j 
%     end
%     dist_ref=10^10;
%     for i=1:length(bg_mtrx_inside)
%        dist=distance(obs_inside(j,3),obs_inside(j,2), ...
%            bg_mtrx_inside(i,2),bg_mtrx_inside(i,1), wgs84km);
%        if dist < dist_ref
%            index=i;
%            dist_ref=dist;
%        end
%     end
%     mu(j)=bg_mtrx(index,3);
% end

disp('Step 2')


%% SOLVING THE SYSTEM TO CALCULATE THE PARAMETERS 

x0 = [10/365,0.1] ;     % initial guess
fun = @(x) root2D (x,mu,C,t2,t1,H,Delta,N_eqk_inside) ;
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter') ;
x = fsolve(fun,x0,options) ;

bgd = x(1) ; % average daily background seismicity rate

v = bgd*365.25 ;  % annual number of background events (events/year)
A = x(2) ;        % estimated simplETAS parameter: earthquake productivity


%%

% The following is a useful test to check if the obtained background (v) is reasonable. 
% It should be roughly comparable to the minimum slope of this curve in some temporal subsets

xcheck = time_inside ;
ycheck = 1:1:N_eqk_inside ;

figure(1)
plot(xcheck,ycheck) ;
