%################################################################
%
% r = scalescape(mimg,maxscale,conv_type)
%
% Create a series of Gaussian blurred images according to a 
% given maximum scale.
%
% Input:    mimg:       image for convolution w/ Gaussian filter
%           maxscale:   maximum scale value
%           conv_type:  convolution type flag
%
% Output:   blurred_imgs:   maxscale blurred images
%
%################################################################

function[blurred_imgs] = scalespace(mimg,maxscale,conv_type)

for scale = 1:1:maxscale
	
	if (scale < 3)
        % Image is unblurred:       
		blurred_imgs(:,:,scale) = mimg;
         
	else
	    % Set values for generating Gaussian filter at given scale:
        [stdd,sizz] = setvalues(scale);
	    kern        = d2gauss(sizz,stdd,1,1,0,1/(stdd*sqrt(2*pi)));	

        c1mimg  = convolve_2(mimg,kern,conv_type);
	    cres    = convolve_2(c1mimg,rot90(kern,3),conv_type);
    
        blurred_imgs(:,:,scale) = cres;
	
	end;
	
end;
	
return;

% Results differ to cantata output for scale >= 3 in the 4th d.p.