function [ lsjoint, lslogz, final_labels, on_nodes, likelihoods ] = trace_Back(DP_table, aux_table, z, node_idx, obsLike)
% It traces back the dinamic programming results assigning labels to each
% position (On, Off). It also computes the joint probability and the log
% evidenct for each line segment. Additionally computes the sum of the
% marginals of all locations along the detected segments
%
% Input
%   - DP_Table: Table that contains the partial results of the Dynamic
%   programing solution
%   - aux_table: table that contains the labels and it will be used to
%   unbound the results
%   - z: Normalization constant
%   - gamma: Marginals at each position along the line
%   - node_idx: Locations along the line to be considered
%
% Output:
%   - lsjoint: Joint probabilty at segment level
%   - lslogz: Evidence of line segments
%   - lsgamma: Sum of the marginals of the positions along segments
%   - final_labels: labels associated with the positions along the line
%   - on_nodes: list of indeces associated with a segment

final_labels = zeros(1, length(node_idx));
[~, final_labels(1,end)] = max(DP_table(:,end));
on_nodes = [];
lsjoint = [];
lslogz = [];
likelihoods = [];
prev = 2; %control start and end line segments -- E.A.
logz_tmp = 0;
ctrl = 0; % used to control end of the line cases
lsId = 1;
for i = length(node_idx)-1:-1:1
   ind = final_labels(i+1);
   final_labels(i) = aux_table(ind, i+1);
   %sum the log of the partial evidences.
   if ctrl == 1 && ind == 1 
       logz_tmp = logz_tmp + log(z(i+1));
%       lsgamma(:,lsId) = lsgamma(:,lsId) + gamma(:,i+1);
       likelihoods(lsId) = likelihoods(lsId) + obsLike(1,i+1);
       on_nodes = [on_nodes; node_idx(i+1)];
   end
   if ind == 1 && prev == 2 %begining of a segment
       endProb = DP_table(ind,i+1);
       logz_tmp = log(z(i+1));
%       lsgamma(:,lsId) = gamma(:,i+1); 
       likelihoods(lsId) = obsLike(1,i+1);
       on_nodes = [on_nodes; node_idx(i+1)];
       ctrl = 1;
   elseif ind == 2 && prev == 1
       startProb = DP_table(ind, i+1);
       lsjoint = [lsjoint; endProb - startProb];
       lslogz = [lslogz; logz_tmp];
       lsId = lsId + 1;
       logz_tmp = 0;
       ctrl = 0;
   elseif ind == 1 && i == 1
       lsjoint = [lsjoint; endProb];
       lslogz = [lslogz; logz_tmp];

   end
   prev = ind;
   %%
end
        
%% E.A. Posterior prob.
lsjoint = flipud(lsjoint); %E.A. To match the line segment output
lslogz = flipud(lslogz); %E.A. To match the line segment output

likelihoods = fliplr(likelihoods);

end

