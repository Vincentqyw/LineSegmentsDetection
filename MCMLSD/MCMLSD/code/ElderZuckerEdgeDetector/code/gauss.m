% gauss.m: Calculates Gaussian pdf values.

function[y] = gauss(x,std)

y = exp(-x.^2/(2*std^2)) / (std*sqrt(2*pi));

return;
