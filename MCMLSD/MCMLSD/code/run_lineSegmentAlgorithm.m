function [out, P] = run_lineSegmentAlgorithm(kernels, kernels_flip, kernel_params, sig_bound, r_res, th_res, img,scale_factor,pp)                                  



frac = 0.05;
%the maximum number of lines that could be extract from the algorithm
maxLine=5000;
LIKE_MODEL = load(['parameters', filesep, 'likelihood_model.mat']);
[P, edgeStruct] = determine_hough_lines_kernel_vote_remove(img, pp, kernels, kernels_flip, maxLine, kernel_params, r_res, th_res, sig_bound, frac,scale_factor);
[begin_points, end_points,ls_likelihoods] = get_all_segments_assoc_edgeremoval( P, edgeStruct, pp, LIKE_MODEL, scale_factor);
ttlLS = size(begin_points,1);
out = zeros(ttlLS,5);  
for dd = 1:ttlLS 
    out(dd,:) = [begin_points(dd,1), begin_points(dd,2), end_points(dd,1), end_points(dd,2),ls_likelihoods(dd)]; %E.A.         
end




