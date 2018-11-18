%###########################################################################
%
% g1scale(g1mag1,g1dir1,g1mag2,g1dir2,g1scale1,scale,noise,b_est)
%
% Purpose: Augments multi-scale Gaussian Gradient maps with significant 
%		estimates at a new scale. Pixels for which the magnitude 
%		of the Gaussian gradient is under threshold in the multi-scale 
%		map (gradient magnitude input 1) but over threshold 
% 		at the new scale (gradient magnitude input 2) are updated 
%		with the gradient magnitude, direction and scale value of 
%		the new scale.
%
% Input:  g1mag1    - Multi-scale Gaussian gradient magnitude image
%         g1dir1    - Multi-scale Gaussian gradient direction image
%         g1mag2    - Gaussian gradient magnitude image at new scale
%         g1dir2    - Gaussian gradient direction image at new scale
%         g1scale1  - Multi-scale scale map
%         g1scale2  - Scale of new gradient estimates
%         noise     - Estimated sensor noise
%         b_est     - Derivatives near boundary estimated by reflecting
%                           intensity function.
% Output: g1mag1    - Integrated multi-scale gradient magnitude map
%         g1dir1    - Integrated multi-scale gradient direction map
%         g1scale1  - Integrated multi-scale scale map
%
%###########################################################################

function[g1mag1,g1dir1,g1sc1] = g1scale(g1mag1,g1dir1,g1mag2,g1dir2,g1sc1,...
                                        scale,noise,b_est)


norms12	= [0.765, 0.199, 0.0499, 0.0125, 0.00312, 0.00078];
thresh 	= 5.6 * noise * norms12(scale);

if ((scale<3) | b_est)
    krad	= 1;
else
    krad 	= ceil(4.6*sqrt(pow2(2*(scale-2))-1));
end
 

if (scale == 1)

    g1mag1 	= zeros(size(g1mag2));
    g1dir1 	= 4 * ones(size(g1mag2));
    g1sc1 	= zeros(size(g1mag2));

    f 			= find(g1mag2 >= thresh);
    g1mag1(f) 	= g1mag2(f);
    g1dir1(f) 	= g1dir2(f);
    g1sc1(f) 	= scale;

else
    sz = size(g1mag2);

	[i,j] 	= meshgrid(krad+1:1:sz(1)-krad,krad+1:1:sz(2)-krad);
    itrans  = i';
    jtrans  = j';
    if ~(isempty(itrans) | isempty(jtrans))
    	K		= sub2ind(size(g1mag2),itrans(:),jtrans(:));

	    smat	= g1sc1(krad+1:end-krad,krad+1:end-krad);
	    mmat	= g1mag2(krad+1:end-krad,krad+1:end-krad);

  	    f 		= find((smat==0) & (mmat>=thresh));

        g1mag1(K(f)) 	= g1mag2(K(f));
        g1dir1(K(f)) 	= g1dir2(K(f));
        g1sc1(K(f)) 	= scale;
    end;
end;

return;
