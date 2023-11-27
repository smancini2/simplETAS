%% From Mancini and Marzocchi (2023), SRL

%% Version of Aug 28, 2023
%
% Acknowledgment: some functions and the main structure of this code were originally written as 
% part of an R code by Dr Stefanie Seif and used for the following paper:

% Seif, S., A. Mignan, J. D. Zechar, M. J. Werner, and S. Wiemer (2017). 
% Estimating ETAS: The effects of truncation, missing data, and model assumptions, 
% J. Geophys. Res. Solid Earth 121, 449-469. https://doi.org/10.1002/2016JB012809.  

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
% THIS FUNCTION SIMULATES simplETAS AFTERSHOCKS %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fin_gen = simulate_aftershocks(parents_init,theta,mag_large,mbin,tmin_aft,tmax_aft)
          
disp('****************************************') ;
disp('Check simplETAS parameters and input file:') ;
fprintf('theta: A=%g, c=%g, alpha=%g, p=%g, gamma=%g, D=%g, q=%g, beta=%g \n',theta(1),theta(2),theta(3),theta(4),theta(5),theta(6),theta(7),theta(8)) ; 

% calculate expected branching ratio - equations (4) from Seif et al. (2017)

if (theta(3) == theta(8))
   br = (theta(1)*theta(8)*(mag_large-theta(9)))/(1-exp(-theta(8)*(mag_large-theta(9)))) ;
else
   br = (theta(1)*theta(8)/(theta(8)-theta(3)))*(1-exp(-(theta(8)-theta(3))*(mag_large-theta(9))))/(1-exp(-theta(8)*(mag_large-theta(9)))) ;
end 
fprintf('\n Expected branching ratio = %g \n', br) ;

% some preliminary checks

if br > 1 && (theta(3) >= theta(8))
   error('Error. The process is explosive! \n') ;
   return
else
    fprintf('\n The process is not explosive \n') ;
end

if isempty(find(isnan(parents_init)))==0
   error('Error. The catalog contains NANs!! \n') ;
   return
else
    fprintf('\n The catalog looks ok \n') ;
end

disp('****************************************') ;

const = theta(1)*exp(theta(3)*(parents_init(:,4)-(theta(9)-mbin/2))) ;
next_gen = parents_init ;
fin_gen = next_gen ;
Naft = poissrnd(const) ;
Nnext = sum(Naft) ;

fprintf('\n Number of aftershocks: %g \n', Nnext) ;

K = 0 ;

while (Nnext > 0.1)
      K = K+1 ;
      disp('****************************************') ;
      fprintf('\n Generation: %g \n', K) ;
      
      taft = repelem(next_gen(:,1), Naft, 1) ;
      index_p = repelem(next_gen(:,7), Naft, 1) ;                       % parent of aftershock
      index_nr = (max(next_gen(:,7))+1:1:max(next_gen(:,7))+Nnext)' ;   % new index for aftershock
      taft = taft + (rand(Nnext,1).^(1/(1-theta(4)))-1) * theta(2) ;    % Omori law
      maft = GRTruncSimul(theta(8),mag_large,theta(9),mbin, Nnext) ;    % truncated G-R law at mag_large
      maft = round(maft,2) ;
      %fprintf('Magnitudes: %g \n', maft) ;
      
      dd = repelem(theta(6)*exp(theta(5)*(next_gen(:,4)-theta(9))), Naft, 1) ;
      RR = sqrt(abs(dd.*(rand(Nnext,1).^(1/(1-theta(7)))-1))) ;         % in km (TO DO: SET A DISTANCE LIMIT INSIDE TESTING POLYGON)
      angle = rand(Nnext,1)*pi*2 ;
      xaft = repelem(next_gen(:,2), Naft, 1) ;
      yaft = repelem(next_gen(:,3), Naft, 1) ;
      %xaft = xaft+(RR*cos(angle))* Km2Lon(yaft) ;   % lon: transform km -> degree
      %yaft = yaft+(RR*sin(angle))* Km2Lat ;         % lat: transform km -> degree
      xaft = xaft+(RR.*cos(angle)) ;                 % location x
      yaft = yaft+(RR.*sin(angle)) ;                 % location y
      
      next_gen = [taft xaft yaft maft ones(size(taft,1),1) K*ones(size(taft,1),1) index_nr index_p] ;

      % for higher order aftershock generations, only keep events which are inside our time window
      [row ~] = find(next_gen(:,1)<tmax_aft & next_gen(:,1)>=tmin_aft) ;
      next_gen = next_gen(row,:) ; clear row ;
      Nnext = size(next_gen,1) ;

      const = theta(1)*exp(theta(3)*(next_gen(:,4)-(theta(9)-mbin/2))) ;
      Naft = poissrnd(const) ;
      Nnext = sum(Naft) ;
      fprintf('\n Number of aftershocks: %g \n', Nnext) ;

      fin_gen = cat(1,fin_gen,next_gen) ;
end

end
