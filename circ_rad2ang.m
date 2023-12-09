function alpha = circ_rad2ang(alpha)

% alpha = circ-rad2ang(alpha)
%   converts values in radians to values in degree
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

alpha = mod(alpha / pi *180,360);

% Modified by adding modulo 360; 
% Added by Nele Schuff (09-08-2023)