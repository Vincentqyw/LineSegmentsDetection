function [final_labels,  begin_points, end_points, begin_point, end_point, posterior, lslength, logratioposterior, ls_probDiff, likelihoods, edgemap] = get_line_segment_DP_edgeremoval( edgeStruct, rho, theta, pp, probStruct, scale_factor)
% Description: Line segment detector that removes the edges associated with
% each detected line segment.

final_labels = [];
begin_points = []; 
end_points = [];
begin_point = []; 
end_point = [];
posterior = [];
lslength = [];
logratioposterior = [];
ls_probDiff = [];

edgemap = edgeStruct.edge;
edge_gradients = edgeStruct.g1dir;
exist_edge_on = probStruct.exist_edge_on;
exist_edge_off = probStruct.exist_edge_off;
ang_dev_on = probStruct.ang_dev_on_hist;
ang_dev_off = probStruct.ang_dev_off_hist;
nexist_edge_on = probStruct.nexist_edge_on;
nexist_edge_off = probStruct.nexist_edge_off;

%the transition probability for the line segments
%it is generated from sampling in the YorkUDB dataset
prob_leave_on = 0.0051/scale_factor;
prob_stay_on =1-prob_leave_on;
prob_leave_off = 0.0014/scale_factor;
prob_stay_off = 1-prob_leave_off;
prob_on = 0.247; 
prob_off = 0.753;
[m, n] = size(edgeStruct.edge);
 
% GET EQUATION OF THE LINE
X_1 = rho * cos(theta) - 2000*cos(theta+pi/2) + pp(1);
X_2 = rho * cos(theta) + 2000*cos(theta+pi/2) + pp(1);
Y_1 = rho * sin(theta) - 2000*sin(theta+pi/2) + pp(2);
Y_2 = rho * sin(theta) + 2000*sin(theta+pi/2) + pp(2);
M = (Y_2-Y_1)/(X_2-X_1 + 1e-7);
b = Y_2 - X_2*M;

% GET BOUNDARIES OF THE LINE IN IMAGE FRAME
p1 = [1, M + b];
p2 = [n, M*n + b];
p3 = [(1-b)/M,1];
p4 = [(m-b)/M, m];

points = [p1;p2;p3;p4];
X = points(:,1);
Y = points(:,2);
        
X_valid = (X <= n).*(X >= 1);
Y_valid = (Y <= m).*(Y >= 1);
p_valid = X_valid.*Y_valid;
points_valid = points(p_valid>0, :); %changed find to logical index
X_l = points_valid(:,1)';
Y_l = points_valid(:,2)';
if isempty(X_l) || isempty(Y_l)
    return;
end
%get the line points
indexXY=sampleLine(rho,theta,n,m,pp);
yy=ceil(indexXY(:,2)+pp(2));
xx=ceil(indexXY(:,1)+pp(1));
index = sub2ind([m,n],yy,xx);

p1 = [X_l(1), Y_l(1)];
p2 = [X_l(2), Y_l(2)];
% FIND POINTS CONSIDERED
node_idx=index;
y_h = edgemap(node_idx);
% y_h(y_h > 0) = 1;
% FIND ORIENTATION DIF AND DIST:
TH = pi-theta;
if TH<0, TH=TH+pi;
elseif TH>pi, TH=TH-pi;
end

edge_gradientsmap=edge_gradients(node_idx);
edgeind=edge_gradientsmap < 0;
edge_gradientsmap(edgeind) = edge_gradientsmap(edgeind) + pi;
edgeind=edge_gradientsmap>pi;
edge_gradientsmap(edgeind)= edge_gradientsmap(edgeind)-pi;

edge_img_locations = [edgeStruct.xzero(node_idx),edgeStruct.yzero(node_idx)];

ang_dev = round((edge_gradientsmap - TH)*180/pi);
ang_dev(ang_dev < -90) = ang_dev(ang_dev < -90) + 180;
ang_dev(ang_dev > 90) = ang_dev(ang_dev > 90) - 180;
range_ang = -90:1:90;
range_dis = 0:0.05:5;
loc_points = [xx, yy];
loc_points(y_h > 0,:) = edge_img_locations(y_h > 0,:);
if size(loc_points,1) == 0
    loc_points = [-1000, -1000];
end


%Get the projection of points along the line defined by p1 and p2
[points_proj, distances] = project_point_to_line_segment_vec(p1,p2,loc_points);
        
distances = interp1(range_dis, range_dis, abs(distances), 'nearest');
range_dis = 0:0.05:1.95;
% Eliminate GROSS OUTLIERS
node_idx(distances > 2) = [];
ang_dev(distances > 2) = [];
y_h(distances > 2) = [];
points_proj(distances > 2,:) = [];
distances(distances > 2) = [];
% SORT:
if abs(M) > 1
    [~, idx] = sort(points_proj(:,2));
else
    [~, idx] = sort(points_proj(:,1));
end
points_proj = points_proj(idx,:);
node_idx = node_idx(idx);
distances = distances(idx);
ang_dev = ang_dev(idx);
y_h = y_h(idx);

ang_dev(y_h == 0) = NaN;
%INIT
begin_point = points_proj(1,:);
end_point = points_proj(end,:);
Psi_t = [prob_stay_on prob_leave_off; prob_leave_on prob_stay_off];
pi_var = [prob_on; prob_off];
[~, obslik,~,DP_table, aux_table, z]=hmmParameters(exist_edge_on,exist_edge_off,ang_dev_on,...
    ang_dev_off,node_idx,Psi_t,pi_var,y_h,distances,range_dis,nexist_edge_on,...
    nexist_edge_off,prob_on,prob_off,prob_stay_on,prob_leave_off,prob_leave_on,...
    prob_stay_off,ang_dev,range_ang);

[lsjoint, lslogz,final_labels, node_on, likelihoods] = trace_Back(DP_table, aux_table, z, node_idx, obslik);
edgemap(node_on) = 0;

final_labels(final_labels == 2) = 0;
% RETURN LINE SEGMENTS
[begin_points,end_points,lslength,startpositions,endpositions]=returnLines(points_proj,final_labels);

%compute p(x=0|y) for each line segment
logjointoff = zeros(size(begin_points,1), 1);
for j = 1:size(begin_points,1) %number of line segments
   for i = startpositions(j):endpositions(j)
       if startpositions(j) == 1
           transition = prob_off;
           past = 0;
       elseif i == startpositions(j)
           transition = prob_leave_on;
           past = DP_table(1,i-1);
       else
           transition = prob_stay_off;
           past = DP_table(2,i-1);
       end
       if y_h(i)
           [~, nn] = min(abs(distances(i) - range_dis));
           [~, mm] = min(abs(ang_dev(i) - range_ang));
           p_y_off_1 = exist_edge_off(nn);
           p_t_off_1 = ang_dev_off(mm);
           loc_off = p_y_off_1*p_t_off_1; %likelihood of y given x = off
           logjointoff(j) = logjointoff(j) + log(loc_off) + log(transition) + past;
       else
           [~, nn] = min(abs(distances(i) - range_dis));
           loc_off = nexist_edge_off(nn);
           logjointoff(j) = logjointoff(j) + log(loc_off) + log(transition) + past;
       end

   end
end

for i = 1:size(logjointoff,1)
    if logjointoff(i) > lslogz(i) || logjointoff(i) > lsjoint(i)
        display('Error');
    end
end
if size(lsjoint,1) > 0
    logratioposterior = lsjoint - logjointoff; 
else
    logratioposterior = 0;
end
end


