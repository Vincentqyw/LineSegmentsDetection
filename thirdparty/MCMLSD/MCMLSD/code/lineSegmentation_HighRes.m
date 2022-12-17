function [lines, fullLines] =lineSegmentation_HighRes(img,kernels, kernels_flip, kernel_params)
%use 640x480 image as a reference. the diagonal size of this image is 800
%if image size 640x480 scale factor = 1
%if image size 1280x960 scale factor = 2 ... etc.
%the lines variable contains the detected line segmentations it arranged as
%[x1 y1 x2 y2 probability]
%The fullLines are the detected lines. It is arranged as [rho theta probability]
pp=[307.551305282635,251.454244960136]; %principle point
sig_bound = 3; %threshold for the line detection
r_res = 0.2; %the theshold for the rho
th_res = 0.002; %the threshold for theta
[m,n,~]=size(img);
scale_factor = sqrt(m^2+n^2)/800; %800 is the diagnal length of 640x480 image
%the kernel of the algorithm is calibrated on 640x480 image
%rescale the principle point
pp(1)=pp(1)*n/640;
pp(2)=pp(2)*m/480;
[lines, fullLines] = run_lineSegmentAlgorithm(kernels, kernels_flip, kernel_params, sig_bound, r_res, th_res, img, scale_factor,pp);                    

end