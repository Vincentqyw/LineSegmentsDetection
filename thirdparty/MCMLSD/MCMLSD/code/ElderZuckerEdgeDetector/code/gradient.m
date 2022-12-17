%############################################################
%
% gradient(maxscale,noise,gauss_a,conv_type,filtpath)
% Computes non-zero gradient in the luminance function by 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[g1mag,g1dir,g1sc] = gradient(maxscale,noise,gauss_a,...
								        conv_type)

fm2 	= '.ascii';

%############
% Initialise:
%############
g1mag	= 0;
g1dir 	= 0;
g1sc	= 0;

for scale = 1:1:maxscale

    % set scale values for filters 
    if (scale == 1)
        g1scaleval = '05';
    else
        g1scaleval = '1';
    end;
    
    mimg	= gauss_a(:,:,scale);

	%#####################################################
    % Compute response of first derivative Gaussian filter 
	%	to the blurred image in an arbitrary direction:
	%#####################################################
    %kern1 	= load(strcat(filtpath,'/gy',g1scaleval,fm2));
    kern1 	= load(strcat('gy',g1scaleval,fm2));
    rc1 	= convolve_2(mimg,kern1,conv_type);
    %kern2 	= load(strcat(filtpath,'/g1x',g1scaleval,fm2));
    kern2 	= load(strcat('g1x',g1scaleval,fm2));
    rc2 	= convolve_2(rc1,kern2,conv_type); % x basis input
    
    %kern3 	= load(strcat(filtpath,'/gx',g1scaleval,fm2));
    kern3 	= load(strcat('gx',g1scaleval,fm2));
    rc3 	= convolve_2(mimg,kern3,conv_type);
    %kern4 	= load(strcat(filtpath,'/g1y',g1scaleval,fm2));
    kern4 	= load(strcat('g1y',g1scaleval,fm2));
    rc4 	= convolve_2(rc3,kern4,conv_type); % y basis input
      
	%###################################################
    % Calculate magnitude and direction of the gradient:
	%###################################################
    [m2,d2] = g1steer(rc2,rc4);
    
    omag = g1mag;
    odir = g1dir;
    osc = g1sc;
    
	%############################################
    % Augment multi-scale Gaussian Gradient maps:
	%############################################
    [g1mag,g1dir,g1sc] = g1scale(g1mag,g1dir,m2,d2,g1sc,scale,noise,0); 
        
end;

return;
