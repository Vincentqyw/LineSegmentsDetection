%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% g2scale - Augments multi-scale Gaussian directional 2nd derivative maps
% with significant estimates at a new scale. Pixels for which the magnitude
% of the 2nd derivative is under threshold in the multi-scale map 
% (2nd derivative input 1) but over threshold at the new scale (2nd derivative
% input 2) are updated with the 2nd derivative magnitude and scale value of the
% new scale.
%
% Input:  g2mag1   - Multi-scale Gaussian directional 2nd derivative image
%         g2mag2   - Gaussian directional 2nd derivative image at new scale
%         g1scale1 - Multi-scale scale map
%         noise    - Estimated sensor noise
%         b_est    - Estimate derivatives near boundaries?
%
% Output: g2mag   - Integrated multi-scale directional 2nd derivative map
%         g2scale - Integrated multi-scale 2nd derivative scale map
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[g2mag1,g2sc1] = g2scale(g2mag1,g2mag2,g2sc1,scale,noise,b_est)

norms12 = [1.873 0.2443 0.0306 0.003817 0.00047715 0.0000596 0.000007455];
thresh 	= 5.2 * noise * norms12(scale);

if ((scale<3) | b_est)
    krad = 1;
else
    krad = ceil(4.6*sqrt(pow2(2*(scale-2))-1.0));
end;

if (scale == 1)

	g2mag1	= zeros(size(g2mag2));
	g2sc1	= zeros(size(g2mag2));

    f 			= find(abs(g2mag2) >= thresh);
    g2sc1(f) 	= scale;
    g2mag1(f) 	= g2mag2(f);

else

   	sz = size(g2mag1); 
	
	[i,j] 	= meshgrid(krad+1:sz(1)-krad,krad+1:sz(2)-krad);
    itrans  = i';
    jtrans  = j';
    if ~(isempty(itrans) | isempty(jtrans))

		K		= sub2ind(size(g2mag1),itrans(:),jtrans(:));
	
       	magmat1 = g2mag1(krad+1:end-krad,krad+1:end-krad);
       	magmat2 = g2mag2(krad+1:end-krad,krad+1:end-krad);
       	f 		= find(abs(magmat1) == 0 & abs(magmat2) >= thresh);
	
       	g2mag1(K(f)) = g2mag2(K(f));
       	g2sc1(K(f))  = scale;
    end;
end;

return;
