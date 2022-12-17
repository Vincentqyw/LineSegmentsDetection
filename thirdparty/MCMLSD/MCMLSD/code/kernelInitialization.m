function [kernels, kernels_flip, kernel_params] =kernelInitialization(img)
[m,n,~]=size(img);
scale_factor = sqrt(m^2+n^2)/800; %800 is the diagnal length of 640x480 image
maxr = floor(400*scale_factor);
sigmax = 0.4;
r_res = 0.2;
th_res = 0.002;
sigmatheta = 5.4*(pi/180);
[kernels, kernels_flip, kernel_params] = precompute_kernels_sparse_res2(maxr, sigmax, sigmatheta, r_res, th_res);
end
