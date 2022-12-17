function [alpha_t, obslik,T,DP_table, aux_table, z]=hmmParameters(exist_edge_on,exist_edge_off,ang_dev_on,...
    ang_dev_off,node_idx,Psi_t,pi_var,y_h,distances,range_dis,nexist_edge_on,...
    nexist_edge_off,prob_on,prob_off,prob_stay_on,prob_leave_off,prob_leave_on,...
    prob_stay_off,ang_dev,range_ang)
node_idx_length=length(node_idx);
DP_table = zeros(2, node_idx_length);
aux_table = zeros(2, node_idx_length);
% aux_table(1,1) = -Inf;
% aux_table(2,1) = -Inf;
z = zeros(node_idx_length,1); % E.A. List of log evidences: sum(log(z)) from 1:t.
alpha_t = zeros(2,node_idx_length); %E.A. 2x1 vector (on/off) of posterior prob.: e.g. alpha_1 = [p(x=on|x); p(x=off|x)]
obslik = zeros(2,node_idx_length); %E.A. 2x1 observation likelihood
if y_h(1)
    [val, nn] = min(abs(distances(1) - range_dis));
    [val, mm] = min(abs(ang_dev(1) - range_ang));
    p_y_on_1 = exist_edge_on(nn);
    p_y_off_1 = exist_edge_off(nn);
    p_t_on_1 = ang_dev_on(mm);
    p_t_off_1 = ang_dev_off(mm);

    loc_on = p_y_on_1*p_t_on_1; %likelihood of y given x = on
    loc_off = p_y_off_1*p_t_off_1; %likelihood of y given x = off

    DP_table(1,1) = log(loc_on) + log(prob_on);
    DP_table(2,1) = log(loc_off) + log(prob_off);

    %Compute prob. of the evidence. E.A.
    obslik(:,1) = [loc_on; loc_off];
    z(1) = sum(obslik(:,1) .* pi_var);
    alpha_t(:,1) = [ loc_on *  pi_var(1); loc_off *  pi_var(2)] / z(1);

%     if (round(sum(alpha_t(:,1))) ~= 1)
%         display('Error');
%     end
            
%% Compute log likelihood ration and evidenct -- E.A             
%             log_like_ratio_list(1) = log((p_y_on_1*p_t_on_1)/(p_y_off_1*p_t_off_1)); %E.A. log likelihood ratio (1st position, edge)
%             evidence(1) = log(p_y_on_1*p_t_on_1*prob_on + p_y_off_1*p_t_off_1*prob_off);
%%      
else
    [val, nn] = min(abs(distances(1) - range_dis));

    loc_on = nexist_edge_on(nn);
    loc_off = nexist_edge_off(nn);

    DP_table(1,1) = log(loc_on) + log(prob_on);
    DP_table(2,1) = log(loc_off) + log(prob_off);


    %Compute prob. of the evidence. E.A.
    %psi = [loc_on; loc_off];
    obslik(:,1) = [loc_on; loc_off];
    z(1) = sum(obslik(:,1) .* pi_var);
    alpha_t(:,1) = [ loc_on *  pi_var(1); loc_off *  pi_var(2)] / z(1);
    %obslik(:,1) = [loc_on; loc_off];

%     if (round(sum(alpha_t(:,1))) ~= 1)
%         display('Error');
%     end           
%             
%% Compute log likelihood ration and evidenct -- E.A             
%             log_like_ratio_list(1) = log(p_y_on_1/p_y_off_1); %E.A. log likelihood ratio (1st position, no edge)  
%             evidence(1) = log(p_y_on_1*prob_on + p_y_off_1*prob_off);
%%      
end

% LOOP THROUGH
for i =2:node_idx_length
    if y_h(i)
        [val, nn] = min(abs(distances(i) - range_dis));
        [val, mm] = min(abs(ang_dev(i) - range_ang));
        p_y_on_1 = exist_edge_on(nn);
        p_y_off_1 = exist_edge_off(nn);
        p_t_on_1 = ang_dev_on(mm);
        p_t_off_1 = ang_dev_off(mm);

        loc_on = p_y_on_1*p_t_on_1; %likelihood of y given x = on
        loc_off = p_y_off_1*p_t_off_1; %likelihood of y given x = off

        [DP_table(1,i), aux_table(1,i)] = max([log(loc_on)  + log(prob_stay_on) + DP_table(1,i-1), log(loc_on) + log(prob_leave_off) + DP_table(2,i-1)]);
        [DP_table(2,i), aux_table(2,i)] = max([log(loc_off) + log(prob_leave_on) + DP_table(1,i-1), log(loc_off) + log(prob_stay_off) + DP_table(2,i-1)]);


        %Compute prob. of the evidence. E.A.            
        %psi = [loc_on; loc_off];
        obslik(:,i) = [loc_on; loc_off];
        predict = Psi_t * alpha_t(:,i-1);
        z(i) = sum(obslik(:,i) .*  predict);
        alpha_t(:,i) = [ loc_on * predict(1); loc_off * predict(2)] / z(i);
        %obslik(:,i) = [loc_on; loc_off];



%         if (round(sum(alpha_t(:,1))) ~= 1)
%             display('Error');
%         end                

%% Compute log likelihood ration and evidenct -- E.A                 
%                 log_like_ratio_list(i) = log((p_y_on_1*p_t_on_1)/(p_y_off_1*p_t_off_1)); %E.A. store likelihood              
%                 evidence(i) = evidence(i-1) + log(p_y_on_1*p_t_on_1*prob_stay_on + p_y_on_1*p_t_on_1*prob_stay_off + ...
%                     p_y_off_1*p_t_off_1*prob_leave_on +  p_y_off_1*p_t_off_1*prob_stay_off); %E.A. log evidence 1..i
%%
   else
        [val, nn] = min(abs(distances(i) - range_dis));

        loc_on = nexist_edge_on(nn);
        loc_off = nexist_edge_off(nn);

        [DP_table(1,i), aux_table(1,i)] = max([log(loc_on)  + log(prob_stay_on) + DP_table(1,i-1), log(loc_on) + log(prob_leave_off) + DP_table(2,i-1)]);
        [DP_table(2,i), aux_table(2,i)] = max([log(loc_off) + log(prob_leave_on) + DP_table(1,i-1), log(loc_off) + log(prob_stay_off) + DP_table(2,i-1)]);

        %Compute prob. of the evidence. E.A.
        obslik(:,i) = [loc_on; loc_off];
        predict = Psi_t * alpha_t(:,i-1);
        z(i) = sum(obslik(:,i) .*  predict);
        alpha_t(:,i) = [ loc_on * predict(1); loc_off * predict(2)] / z(i);

%         if (round(sum(alpha_t(:,1))) ~= 1)
%             display('Error');
%         end

%% Compute log likelihood ration and evidence --  E.A                 
%                 log_like_ratio_list(i) = log(p_y_on_1/p_y_off_1); %E.A. store likelihood
%                 evidence(i) = evidence(i-1) + log(p_y_on_1*prob_stay_on + p_y_on_1*prob_leave_off + ...
%                     p_y_off_1*prob_leave_on + p_y_off_1*prob_stay_off); %E.A. log evidence 1..i
%%                
    end
end
T = size(DP_table,2);
end