%% round2frac() Rounds the input data to the nearest fraction
%
%   Description: 
%
%       Takes 'data' a matrix of floats and rounds it to the nearest 'frac'
%
%   Author:
%
%       Ron Tal 
%
%   Date Created:
%
%       November 4, 2008
%
function [round_data] = round2frac(data, frac)
%     if frac > 1, error('Fraction must be smaller than 1'); end
    round_data = data.*(1/frac);
    round_data = round(round_data);
    round_data = round_data.*frac;
end