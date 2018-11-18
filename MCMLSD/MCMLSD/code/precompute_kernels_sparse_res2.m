function [ kernels, kernels_flip, kernel_params] = precompute_kernels_sparse_res2(r_max, SIGMA_X, SIGMA_TH, r_res, th_res)
%PRECOMPUTE_KERNELS Precomputes the kernels for Hough Transform
%
%   Inputs:
%               r_max:      Maximum value r_th can take
%               SIGMA_X:    Variance of the distribution of position
%               SIGMA_TH:   Variance of the distribution of angular
%                           uncertainty
%
%   Output:
%               kernels:    Precomputed kernels for the ranges
%
%   Author:
%               Ron Tal
%
%   Date:
%               January 18, 2009
%
%   Description:
%               This function precomputes the kernels to be used in the
%               Hough Transform method. The kernels are computed for the
%               given range at the given resolutions, using the given sigma
%               values.
%
%   Recommended value for SIGMA_TH: 0.1308
%   Recommended value for SIGMA_X:  0.3933
%
    kernels = cell(1, r_max);
    kernels_flip = cell(1, r_max);
    kernel_params = cell(1, r_max);
    for e = -r_max:1:r_max
        [kern,ind]=computeKernel(e,r_max,SIGMA_X, SIGMA_TH, r_res, th_res);
        kernels(1, ind) = {kern};
        kernels_flip(1, ind) = {kern};
        kernel_params(1, ind) = {max(kern)};
    end  
end
