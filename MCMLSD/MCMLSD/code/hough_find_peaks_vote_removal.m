%% Function P = hough_find_peaks_with_near_peak_supression(map, num_lines, NHOOD)
%           Input
%                   map                    Hough Map
%                   num_lines              Maximum number of lines
%                   NHOOD                  [R, TH] neighborhood to be
%                                          surpressed
%
%           Output
%                   P                      Hough Peaks
%
%
function [P] = hough_find_peaks_vote_removal(map, num, r_res, th_res, edgeStruct, kernels, kernels_flip, pp, kernel_params, sig_th, sig_r, sig_bound, frac, scale_factor)
    
    % Loop until no more peaks or number of desired peaks have been
    % selected
    COUNT = 0;
    MAX_R = 400*scale_factor;
    P =[];
    edgeMap = edgeStruct.edge;
    iter = 1; 
    [m,n]=size(edgeMap);
    while 1
        iter = iter + 1;
       % Find max peak:
        [maxmap,ind]=max(map(:));
        [r,c] = ind2sub(size(map), ind);     
       if length(r) > 1, r = r(1); c = c(1); end
       %here are some stopping criteria 
       if isempty(r)
           return; 
       end
       [rho, theta] = mat2hough_scale(r, c, MAX_R, r_res, th_res);
       if ismember([rho, theta, maxmap], P, 'rows', 'legacy')
           return; 
       end
       if COUNT > 0
           if (maxmap/P(1,3)) < frac 
               return; 
           end
       end
       % Save max peak
       P = [P; [rho theta maxmap]];
       COUNT = COUNT + 1;
        
       % If reached number of lines, return
       if size(P, 1) == num
           return; 
       end
       
       %there following block of code identify the location of edges that
       %need to be removed from the hough map and edge map
        indexXY=sampleLine(rho,theta,n,m,pp);
        yy=ceil(indexXY(:,2)+pp(2));
        xx=ceil(indexXY(:,1)+pp(1));
        edge_index = sub2ind([m,n],yy,xx);
        X_1 = rho * cos(theta) - 2000*cos(theta+pi/2) + pp(1);
        X_2 = rho * cos(theta) + 2000*cos(theta+pi/2) + pp(1);
        Y_1 = rho * sin(theta) - 2000*sin(theta+pi/2) + pp(2);
        Y_2 = rho * sin(theta) + 2000*sin(theta+pi/2) + pp(2);
        M = (Y_2-Y_1)/(X_2-X_1+1e-7);
       b = Y_2 - X_2*M;
       line_eq = [M, -1, b];

       edge_img_locations = [edgeStruct.xzero(edge_index),edgeStruct.yzero(edge_index), ones(length(edge_index(:)),1)];
       edge_gradients = edgeStruct.g1dir(edge_index);
       
       edge_gradients2 = edge_gradients;
       indexless0=edge_gradients2 < 0;
       indexgreaterpi=edge_gradients2>pi;
       edge_gradients2(indexless0) = edge_gradients2(indexless0) + pi;
       edge_gradients2(indexgreaterpi)= edge_gradients2(indexgreaterpi)-pi;
       
       edge_gradients2 = pi - edge_gradients2;
       edgetetalesspi=edge_gradients2 - theta<-pi/2;
       edge_gradients2(edgetetalesspi)= edge_gradients2(edgetetalesspi)+pi;
       edgethetagreaterpi=edge_gradients2-theta>pi/2;
       edge_gradients2(edgethetagreaterpi)= edge_gradients2(edgethetagreaterpi)-pi;
       orientationDiff = abs(edge_gradients2 - theta);
       distances =  (abs(line_eq*edge_img_locations')'/norm(line_eq(1:2)));
       
       distanceBoundary=(distances <= sig_bound*sig_r);
        edge_img_locations_remove  = edge_img_locations(distanceBoundary, 1:2);
        edge_gradients_remove = edge_gradients(distanceBoundary);
              
        orientationDiff = orientationDiff(distanceBoundary);
        edge_index = edge_index(distanceBoundary);
        orientatinDiffboundary=orientationDiff <=sig_bound*sig_th;
        edge_img_locations_remove = edge_img_locations_remove(orientatinDiffboundary, :);
        edge_gradients_remove = edge_gradients_remove(orientatinDiffboundary);
        edge_index = edge_index(orientatinDiffboundary);
       
        if isempty(edge_img_locations_remove) 
            break; 
        end
        edgeMap(edge_index) = 0;
        %remove the detected lines from the hough map
        map=mexRemoveVotes_v3_scale(edge_img_locations_remove, edge_gradients_remove, kernels, kernels_flip, r_res, th_res, pp, kernel_params, map, scale_factor);

        %map = mexRemoveVotes_v3_scale(edge_img_locations_remove, edge_gradients_remove, kernels, kernels_flip, r_res, th_res, pp, kernel_params, map, scale_factor);

    end
end

   