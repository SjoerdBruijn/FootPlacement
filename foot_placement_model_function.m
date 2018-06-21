
function [OUT,intermediates]=foot_placement_model_function(COM,Rfoot,Lfoot,rhs,lhs,fsopto,predict_foot_sample,origin_sample,pred_samples,order,removeorigin,centerdata)
% function to calculate foot placement model and % explained variance.
% returns output structure with variables and suggested figure titles as
% well as y-labels



% include; -intermediates, n# of steps excluded, steptimes
%% actual calculations. Start with setting things up.

% first find some strides that need to be excluded (any that are too long)
st_L = diff(lhs);
st_R = diff(rhs);
ignore_L = find(st_L>(median(st_L)+0.3*median(st_L))|st_L<(median(st_L)-0.3*median(st_L)))';
ignore_R = find(st_R>(median(st_R)+0.3*median(st_R))|st_R<(median(st_R)-0.3*median(st_R)))';

st_L(ignore_L)              = [];
OUT.stride_time.data        = nanmean(st_L./fsopto);
OUT.stride_time_var.data    = nanstd(st_L./fsopto);

OUT.stride_time.titel       ='Stride time';
OUT.stride_time.ylabel      ='Stride time [s]';

OUT.stride_time_var.titel   ='Stride time variability';
OUT.stride_time_var.ylabel  ='Stride time variability [s]';
%%
% calculate velocities
COM_vel = calc_derivative(COM,fsopto);
COM_acc = calc_derivative(COM_vel,fsopto);

% time normalize; left variables
COM_L       = normalizetimebase(COM,lhs);
COM_L_vel   = normalizetimebase(COM_vel,lhs);
COM_L_acc   = normalizetimebase(COM_acc,lhs);
foot_L      = normalizetimebase(Lfoot,lhs);
origin_L    = normalizetimebase(Rfoot,lhs);

% do the above two things also for right
COM_R       = normalizetimebase(COM,rhs);
COM_R_vel   = normalizetimebase(COM_vel,rhs);
COM_R_acc   = normalizetimebase(COM_acc,rhs);
foot_R      = normalizetimebase(Rfoot,rhs);
origin_R    = normalizetimebase(Lfoot,rhs);

% remove last two colums, and any strides that were too long/short.
COM_L(:,[ignore_L end-1:end])       = nan;
COM_L_vel(:,[ignore_L end-1:end])   = nan;
COM_L_acc(:,[ignore_L end-1:end])   = nan;
foot_L(:,[ignore_L end-1:end])      = nan;
origin_L(:,[ignore_L end-1:end])    = nan;
% remove last two colums
COM_R(:,[ignore_R end-1:end])       = nan;
COM_R_vel(:,[ignore_R end-1:end])   = nan;
COM_R_acc(:,[ignore_R end-1:end])   = nan;
foot_R(:,[ignore_R end-1:end])      = nan;
origin_R(:,[ignore_R end-1:end])    = nan;
%% save outcomes
OUT.COM_var.data     = nanstd(COM_L,1,2);
OUT.COM_vel_var.data = nanstd(COM_L_vel,1,2);
OUT.COM_acc_var.data = nanstd(COM_L_acc,1,2);

OUT.COM_var.ylabel      = 'CoM variability [m]';
OUT.COM_vel_var.ylabel  = 'CoM velocity variability [m/s]';
OUT.COM_acc_var.ylabel  = 'CoM acceleration variability [m/s^2]';

OUT.COM_var.titel       = 'CoM variability';
OUT.COM_vel_var.titel   = 'CoM velocity variability';
OUT.COM_acc_var.titel   = 'CoM acceleration variability';

OUT.COM.data     = nanmean(COM_L,2);
OUT.COM_vel.data = nanmean(COM_L_vel,2);
OUT.COM_acc.data = nanmean(COM_L_acc,2);

OUT.COM.ylabel      = 'CoM position [m]';
OUT.COM_vel.ylabel  = 'CoM velocity [m/s]';
OUT.COM_acc.ylabel  = 'CoM acceleration [m/s^2]';

OUT.COM.titel       = 'CoM position';
OUT.COM_vel.titel   = 'CoM velocity';
OUT.COM_acc.titel   = 'CoM acceleration';

%%
foot_L       = foot_L(predict_foot_sample,:)';
origin_L     = origin_L(origin_sample,:)';
foot_R       = foot_R(predict_foot_sample,:)';
origin_R     = origin_R(origin_sample,:)';

if removeorigin
    foot_L       = foot_L -origin_L; % substract origin, which is the other foot
    foot_R       = foot_R -origin_R ;
end
OUT.SW.data         = nanmean(abs([foot_L; foot_R]));
OUT.SW_var.data     = nanstd(abs([foot_L ;foot_R]));

OUT.SW.ylabel       = 'Stepwidth [m]';
OUT.SW_var.ylabel   = 'Stepwidth variability[m]';
OUT.SW.titel        = 'Stepwidth';
OUT.SW_var.titel    = 'Stepwidth variability';

if centerdata
    foot_L       = foot_L -nanmean(foot_L); % substract mean
    foot_R       = foot_R -nanmean(foot_R);
end

%%
L_jac           = zeros(length(pred_samples),order);
R_jac           = zeros(length(pred_samples),order);
combined_jac    = zeros(length(pred_samples),order);
for i_pred_sample = pred_samples
    %  get correct sample for all variables
    COM_L_sample        = COM_L(i_pred_sample,:)';
    COM_L_vel_sample    = COM_L_vel(i_pred_sample,:)';
    COM_L_acc_sample    = COM_L_acc(i_pred_sample,:)';
    
    % get correct sample
    COM_R_sample        = COM_R(i_pred_sample,:)';
    COM_R_vel_sample    = COM_R_vel(i_pred_sample,:)';
    COM_R_acc_sample    = COM_R_acc(i_pred_sample,:)';
    
    % predictors
    pred_Lstance = [COM_L_sample  COM_L_vel_sample  COM_L_acc_sample ];
    pred_Rstance = [COM_R_sample  COM_R_vel_sample  COM_R_acc_sample ];
    
    % remove origin and means.
    if removeorigin
        pred_Lstance(:,1)   = pred_Lstance(:,1)-origin_L;
        pred_Rstance(:,1)   = pred_Rstance(:,1)-origin_R;
        OUT.pred1_Leftstance.data(i_pred_sample) = nanmean(pred_Lstance(:,1)); 
        OUT.pred1_Leftstance.ylabel = 'Left: predictor 1';
        OUT.pred1_Leftstance.titel  = 'Left: predictor 1';
        OUT.pred1_Rightstance.data(i_pred_sample) = nanmean(pred_Lstance(:,1)); 
        OUT.pred1_Rightstance.ylabel    = 'right: predictor 1';
        OUT.pred1_Rightstance.titel     = 'right: predictor 1';
    end
    if centerdata
        pred_Lstance        = pred_Lstance-repmat(nanmean(pred_Lstance),size(pred_Lstance,1),1);
        pred_Rstance        = pred_Rstance-repmat(nanmean(pred_Rstance),size(pred_Rstance,1),1);
    end
    foot_L_sample   = foot_L;
    foot_R_sample   = foot_R;
    % save predictors
    OUT.var_pre1dLeftstance.data(i_pred_sample)    = nanstd(pred_Lstance(:,1));
    OUT.var_pre2dLeftstance.data(i_pred_sample)    = nanstd(pred_Lstance(:,2));
    
    OUT.var_pre1dLeftstance.ylabel   = 'Left: Variability of predictor 1 [m]';
    OUT.var_pre2dLeftstance.ylabel   = 'Left: Variability of predictor 2 [m/s]';
    OUT.var_pre1dLeftstance.titel    = 'Left: Variability of predictor 1';
    OUT.var_pre2dLeftstance.titel    = 'Left: Variability of predictor 2';
    
    OUT.var_pre1dRightstance.data(i_pred_sample)    = nanstd(pred_Rstance(:,1));
    OUT.var_pre2dRightstance.data(i_pred_sample)    = nanstd(pred_Rstance(:,2));
    OUT.var_pre1dRightstance.ylabel   = 'Right: Variability of predictor 1 [m]';
    OUT.var_pre2dRightstance.ylabel   = 'Right: Variability of predictor 2 [m/s]';
    OUT.var_pre1dRightstance.titel    = 'Right: Variability of predictor 1';
    OUT.var_pre2dRightstance.titel    = 'Right: Variability of predictor 2';
    
    %% calculations for right
    tmp                                 = [foot_R_sample pred_Rstance];
    ind_R                               = 1:size(tmp,1);
    foot_R_sample(isnan(sum(tmp,2)))    = [];
    pred_Rstance(isnan(sum(tmp,2)),:)   = [];
    ind_R(isnan(sum(tmp,2)))            = [];
    R_jac(i_pred_sample,:)              = foot_R_sample'/(pred_Rstance(:,1:order)'); % this is the 'jacobian'; as in Wang & srinivasavan
    OUT.Right_pct.data(i_pred_sample)   = corr(foot_R_sample,(R_jac(i_pred_sample,:)*(pred_Rstance(:,1:order)'))','rows','complete')^2; 
    OUT.Right_coeff1.data(i_pred_sample)= R_jac(i_pred_sample,1);
    OUT.Right_coeff2.data(i_pred_sample)= R_jac(i_pred_sample,2);
    OUT.Right_N.data(i_pred_sample)     = size(pred_Rstance,1);
    if i_pred_sample==1
        OUT.Right_foot_var.data      = std(foot_R_sample);
        OUT.Right_foot_var.titel     = 'Right: variance of outcome';
        OUT.Right_foot_var.ylabel    = 'Right: variance of outcome (m)';
    end
    
    OUT.Right_pct.titel     = 'Right: % explained variance';
    OUT.Right_coeff1.titel 	= 'Right: Coefficient 1';
    OUT.Right_coeff2.titel  = 'Right: Coefficient 2';
    OUT.Right_N.titel       = 'Right: number of datapoints included';
    
    OUT.Right_pct.ylabel     = 'Right: % explained variance [%]';
    OUT.Right_coeff1.ylabel  = 'Right: Coefficient 1';
    OUT.Right_coeff2.ylabel  = 'Right: Coefficient 2';
    OUT.Right_N.ylabel       = 'Right: number of datapoints included';
    intermediates.error_right(i_pred_sample,ind_R)   = foot_R_sample-(R_jac(i_pred_sample,:)*(pred_Rstance(:,1:order)'))';
    
    %% for left
    tmp                                 = [foot_L_sample pred_Lstance];
    ind_L                               = 1:size(tmp,1);
    foot_L_sample(isnan(sum(tmp,2)))    = [];
    pred_Lstance(isnan(sum(tmp,2)),:)   = [];
    ind_L(isnan(sum(tmp,2)))            = [];
    L_jac(i_pred_sample,:)              = foot_L_sample'/(pred_Lstance(:,1:order)'); % this is the 'jacobian'; as in Wang & srinivasavan
    OUT.Left_pct.data(i_pred_sample)    = corr(foot_L_sample,(R_jac(i_pred_sample,:)*(pred_Lstance(:,1:order)'))','rows','complete')^2; 
    OUT.Left_coeff1.data(i_pred_sample) = L_jac(i_pred_sample,1);
    OUT.Left_coeff2.data(i_pred_sample) = L_jac(i_pred_sample,2);
    OUT.Left_N.data(i_pred_sample)      = size(pred_Lstance,1);
    if i_pred_sample==1
        OUT.Left_foot_var.data      = std(foot_L_sample);
        OUT.Left_foot_var.titel     = 'Left: variance of outcome';
        OUT.Left_foot_var.ylabel    = 'Left: variance of outcome (m)';
        
    end
    OUT.Left_pct.titel     = 'Left: % explained variance';
    OUT.Left_coeff1.titel  = 'Left: Coefficient 1';
    OUT.Left_coeff2.titel  = 'Left: Coefficient 2';
    OUT.Left_N.titel       = 'Left: number of datapoints included';
    
    OUT.Left_pct.ylabel     = 'Left: % explained variance [%]';
    OUT.Left_coeff1.ylabel  = 'Left: Coefficient 1';
    OUT.Left_coeff2.ylabel  = 'Left: Coefficient 2';
    OUT.Left_N.ylabel       = 'Left: number of datapoints included';
    intermediates.error_left(i_pred_sample,ind_L)   = foot_L_sample-(L_jac(i_pred_sample,:)*(pred_Lstance(:,1:order)'))';
    %% combined
    foot_combined                       = [foot_R_sample ;foot_L_sample];
    pred_combined                       = [pred_Rstance(:,1:order); pred_Lstance(:,1:order)];
    tmp                                 = [foot_combined pred_combined];
    ind_combined                        = 1:size(tmp,1);
    foot_combined(isnan(sum(tmp,2)))    = [];
    pred_combined(isnan(sum(tmp,2)),:)  = [];
    ind_combined(isnan(sum(tmp,2)))     = [];
    combined_jac(i_pred_sample,:)           = foot_combined'/pred_combined'; % this is the 'jacobian'; as in Wang & srinivasavan
    OUT.Combined_pct.data(i_pred_sample)    = corr(foot_combined,(combined_jac(i_pred_sample,:)*pred_combined')','rows','complete')^2;
    OUT.Combined_coeff1.data(i_pred_sample) = combined_jac(i_pred_sample,1);
    OUT.Combined_coeff2.data(i_pred_sample) = combined_jac(i_pred_sample,2);
    OUT.Combined_N.data(i_pred_sample)      = size(foot_combined,1);
    if i_pred_sample==1
        OUT.Combined_foot_var.data      = std(foot_combined);
        OUT.Combined_foot_var.titel     = 'Combined: variance of outcome';
        OUT.Combined_foot_var.ylabel    = 'Combined: variance of outcome (m)';
    end
    OUT.Combined_pct.titel     = 'Combined: % explained variance';
    OUT.Combined_coeff1.titel  = 'Combined: Coefficient 1';
    OUT.Combined_coeff2.titel  = 'Combined: Coefficient 2';
    OUT.Combined_N.titel       = 'Combined: number of datapoints included';
    
    OUT.Combined_pct.ylabel     = 'Combined: % explained variance [%]';
    OUT.Combined_coeff1.ylabel  = 'Combined: Coefficient 1';
    OUT.Combined_coeff2.ylabel  = 'Combined: Coefficient 2';
    OUT.Combined_N.ylabel       = 'Combined: number of datapoints included';
    intermediates.error_combined(i_pred_sample,ind_combined)   = foot_combined-(combined_jac(i_pred_sample,:)*pred_combined')';
end
intermediates.error_right(intermediates.error_right==0)         = nan;
intermediates.error_left(intermediates.error_left==0)           = nan;
intermediates.error_combined(intermediates.error_combined==0)   = nan;

