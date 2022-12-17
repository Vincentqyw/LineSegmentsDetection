%##########################################################################
%
% [stdd,sizz] = setvalues(scale)
%
% Set values for generating 2d Gaussian filter according to 
% input scale.
%
% Input:    scale
% 
% Output:   stdd:       std dev'n along width (along height assumed to be 1)
%           size:       width of output = #columns (height assumed to be 1)
%
%##########################################################################

function[stdd,sizz] = setvalues(scale);

if (scale < 3)
    stdd    = 1;
else
    stdd    = sqrt((2^(scale-2))^2 - 1);
end;

sizz    = 2 * ceil(4.6 * stdd) + 1;

return;