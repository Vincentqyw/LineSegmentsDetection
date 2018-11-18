
%######################################################################
% 
% g1steer - Computes magnitude and direction of the gradient  of the 
% luminance function based on x and y basis functions for 1st Gaussian 
% derivative.
% 
%       
%        Input:  g1x - X basis for Gaussian gradient
%                g1y - Y basis for Gaussian gradient
%        Output: g1mag - Gradient Magnitude Estimate
%                g1dir - Gradient Direction Estimate
%
%######################################################################


function[g1mag,g1dir] = g1steer(g1x,g1y)

g1mag	= zeros(size(g1x));
g1dir	= 4 * ones(size(g1x));

f 	= find((g1x ~= 0) & (g1y ~= 0));

g1dir(f) = atan2(-g1y(f),g1x(f));
g1mag(f) = sqrt(g1y(f).^2 + g1x(f).^2);

return;
