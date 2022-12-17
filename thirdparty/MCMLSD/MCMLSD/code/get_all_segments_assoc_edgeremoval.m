function [ begin_points, end_points,ls_likelihoods] = get_all_segments_assoc_edgeremoval( line_data, edgeStruct, pp, LIKE_MODEL,scale_factor)

begin_points = [];
end_points = [];
ls_likelihoods = [];
for i = 1:size(line_data, 1)
    [~, begin_pointsi, end_pointsi, ~,~, ~, ~, ~, ~, likelihoods, edgeMap] = get_line_segment_DP_edgeremoval( edgeStruct, line_data(i,1), line_data(i,2), pp, LIKE_MODEL, scale_factor); 
    for j = 1:size(likelihoods,2)
            begin_points = [begin_points; begin_pointsi(j,:)];
            end_points = [end_points; end_pointsi(j,:)]; 
            ls_likelihoods = [ls_likelihoods; likelihoods(j)];
    end
    edgeStruct.edge = edgeMap;
end

end

