
function [mag] = GRTruncSimul(beta, Mmax, Mc, M_BIN, nr)

cdf = rand(nr,1) ;
mag = log(-cdf*(-exp(-beta*(Mmax-Mc+M_BIN/2))+1)+1)/(-beta)+ Mc - M_BIN/2 ;
