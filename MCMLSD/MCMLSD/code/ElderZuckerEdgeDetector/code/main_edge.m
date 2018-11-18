%############################################################
%
% main_edge.m:  Main function called from edge GUI interface.
%
%############################################################

function[OutputData] = main_edge(img,scale,noise,edgew,...
                        conss,congrd,subpixelflag)

%#############
% Input Image:
%#############
if (size(img,3)~=1)
    img = rgb2gray(img);
end;
img=double(img);


%########################
% Compose blurred images:
%########################
gauss_imgs  = scalespace(img,scale,conss);

%########################
% Calculate gradient map:
%########################
[g1_mag,g1_dir,g1_sc] = gradient(scale,noise,gauss_imgs,congrd);

%##############################
% Calculate 2nd derivative map:
%##############################
[g2_mag,g2_sc,g2_all] = derivative2nd(g1_dir,scale,noise,gauss_imgs,congrd);


% save input.mat img g1_mag g1_dir g1_sc g2_mag g2_sc g2_all noise edgew subpixelflag scale
%inparTest = checkInParam(img,g1_mag,g1_dir,g1_sc,g2_mag,g2_sc,g2_all,noise,...
%                           edgew,subpixelflag,scale);
%if inparTest > 0
%   display('stop');
%end
%####################
% Calculate edge map:
%####################
clear find_edges
[edge_map,blur_map,dark_map,light_map,xzero1_map,yzero1_map,xzero2_map,yzero2_map] = ...
        find_edges(img,g1_mag,g1_dir,g1_sc,g2_mag,g2_sc,g2_all,noise,...
                            edgew,subpixelflag,scale);
clear find_edges
                        
% outparTest = checkOutParam(edge_map,blur_map,dark_map,light_map,xzero1_map,yzero1_map,xzero2_map,yzero2_map);                        
% if outparTest > 0
%    display('stop');
% end

%save output.mat edge_map blur_map dark_map light_map xzero1_map yzero1_map xzero2_map yzero2_map                         
% Branka's code leaves dark and light as double which need rounding:
dark_map    = round(dark_map);
light_map    = round(light_map);
    

%#######################
% Output data in struct:
%#######################
% OutputData = struct('edge',edge_map,'dark_map',dark_map,...
%                     'light_map',light_map,'blur_map',blur_map,...
%                     'g1mag',g1_mag,'g1dir',g1_dir,'g1scale',g1_sc,...
%                     'g2mag',g2_mag,'g2scale',g2_sc,...
%                     'x1zero',xzero1_map,'y1zero',yzero1_map,...
%                     'x2zero',xzero2_map,'y2zero',yzero2_map,...
%                     'noise',noise);
                
OutputData = struct('edge',edge_map,'blur',blur_map,...
                    'dark',dark_map,'light',light_map,...
                    'g1mag',g1_mag,'g1dir',g1_dir,'g1scale',g1_sc,...
                    'g2mag',g2_mag,'g2scale',g2_sc,'g2_all',g2_all,...
                    'noise',noise,'xzero',xzero1_map,'yzero',yzero1_map);
                
                
                
return;
