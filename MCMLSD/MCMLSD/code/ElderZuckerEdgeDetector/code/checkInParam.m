function [out] = checkInput(img1,g1_mag1,g1_dir1,g1_sc1,g2_mag1,g2_sc1,g2_all1,noise1,...
                            EDGE_W1,subpixelflag1,maxscale1)
 
load('input.mat');
out = 0;
if ~isequal(img1, img)
    display('error img');
    out = out + 1;
end
if ~isequal(g1_mag1, g1_mag)
    display('error g1_mag');
    out = out + 1;
end
if ~isequal(g1_dir1, g1_dir)
    display('error g1_dir');
    out = out + 1;
end
if ~isequal(g1_sc1, g1_sc)
    display('error g1_sc');
    out = out + 1;
end
if ~isequal(g2_mag1, g2_mag)
    display('error g2_mag');
    out = out + 1;
end
if ~isequal(g2_sc1, g2_sc)
    display('error g2_sc');
    out = out + 1;
end
if ~isequal(g2_all1, g2_all)
    display('error g2_all');
    out = out + 1;
end
if noise1 ~= noise
   display('error noise'); 
   out = out + 1;
end
if EDGE_W1 ~= edgew
   display('error EDGE_W'); 
   out = out + 1;
end
if subpixelflag1 ~= subpixelflag
   display('error subpixelflag');
   out = out + 1;
end
if maxscale1 ~= scale
   display('error scale'); 
   out = out + 1;
end
                        
end