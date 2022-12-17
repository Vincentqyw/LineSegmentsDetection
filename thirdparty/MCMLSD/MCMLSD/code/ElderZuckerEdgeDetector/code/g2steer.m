%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% g2steer - Computes the second Gaussian derivative of the luminance function
% in specified direction ( normally the Gradient direction). The three input
% basis function are used to steer the derivative. The units of the direction
% map are radians and the range is between -pi and pi. A value of -4 indicate
% that no direction was measurable. The 2nd derivative is taken only for valid
% directions.
%
% Input:  g2x   - Response map for 1st G2 basis function
%         g2y   - Response map for 2nd G2 basis function
%         g2xy  - Response map for 3rd G2 basis function
%         g1dir - Luminance gradient direction map
%
% Output: g2 - Second derivative response map
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[g2] = g2steer(g2x,g2y,g2xy,g1dir)

g2	    = zeros(size(g2x));

f 		= find((g1dir ~= 4.0) & (g2x ~= 0) & (g2xy ~= 0) & (g2y ~= 0));
cdir 	= cos(2*g1dir);
sdir    = sin(2*g1dir);

g2(f) 	= 0.5 * (1 + cdir(f)) .* g2x(f) - sdir(f) .* g2xy(f) + ...
          0.5 * (1 - cdir(f)) .* g2y(f);

return;
