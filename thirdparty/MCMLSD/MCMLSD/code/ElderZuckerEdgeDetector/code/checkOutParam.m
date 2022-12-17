function  [out] = checkOutParam(edge_map1,blur_map1,dark_map1,light_map1,...
    xzero1_map1,yzero1_map1,xzero2_map1,yzero2_map1) 
                           
load('output.mat');
out = 0;
if ~isequal(edge_map1, edge_map)
    display('error edge_map');
    out = out + 1;
end
if ~isequal(blur_map1, blur_map)
    display('error blur_map');
    out = out + 1;
end
if ~isequal(dark_map1, dark_map)
    display('error dark_map');
    out = out + 1;
end
if ~isequal(light_map1, light_map)
    display('error light_map');
    out = out + 1;
end
if ~isequal(xzero1_map1, xzero1_map)
    display('error xzero1_map');
    out = out + 1;
end
if ~isequal(yzero1_map1, yzero1_map)
    display('error yzero1_map');
    out = out + 1;
end
if ~isequal(xzero2_map1, xzero2_map)
    display('error xzero2_map');
    out = out + 1;
end
if ~isequal(yzero2_map1, yzero2_map)
    display('error yzero2_map');
    out = out + 1;
end      

end