clc
clear
close all
for i = 1:9
    clear mex
   edgeStruct = elderEdge();
   
   figure;
   imshow(edgeStruct.edge_map);
   i
   length(find(edgeStruct.edge_map))  
   diffs{i} = edgeStruct.edge_map;
end