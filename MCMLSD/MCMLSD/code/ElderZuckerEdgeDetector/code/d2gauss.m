%################################################################
%
% d2gauss.m returns a 2-d Gaussian filter with kernal attributes:
%           size:       n1*n2 
%           theta:      CCW-angle tkernat filter rotated
%           sigma1:     standard deviation of 1st gaussian
%           sigma2:     standard deviation of 2nd gaussian
%
%################################################################

function[kern] = d2gauss(n1,std1,n2,std2,theta,max1);

[I,J]   = meshgrid(1:1:n2,1:1:n1);
It      = I - (n2+1)/2;
Jt      = J - (n1+1)/2;

u1      = cos(theta)*Jt' - sin(theta)*It';
u2      = sin(theta)*Jt' + cos(theta)*It';

kern    = gauss(u1,std1) .* gauss(u2,std2);

%############################################
% Normalise the kernal and confine to limits:
%############################################
kern    = kern / sqrt(sum(sum(kern.*kern)));
max2    = max(max(kern));
kern    = kern / (max2/max1);

return;