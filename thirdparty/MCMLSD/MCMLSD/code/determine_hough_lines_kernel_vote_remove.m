function [ lines, edgeStruct ] = determine_hough_lines_kernel_vote_remove( im, pp, kernels, kernels_flip, num, kernel_params, r_res, th_res, sig_bound, frac,scale_factor)

                                                                   
%edge detection
edgeStruct = main_edge(im,5,1,10,1,0,1);

%find the pixel location that has edge value
edge_index = find(edgeStruct.edge);
edge_img_locations = [edgeStruct.xzero(edge_index),edgeStruct.yzero(edge_index)];
edge_gradients = edgeStruct.g1dir(edge_index);
%create a hough map

map = mexVoteEdges_v3_scale(edge_img_locations, edge_gradients, kernels, kernels_flip, r_res,th_res, pp, kernel_params,scale_factor);

%iteratively extract lines from hough map
[P] = hough_find_peaks_vote_removal(map, num, r_res, th_res, edgeStruct, kernels, kernels_flip, pp, kernel_params, 0.0917, 0.3934, sig_bound, frac, scale_factor);
%the lines consist of [rho, theta, strength]
lines = P(:, 1:3);
end
