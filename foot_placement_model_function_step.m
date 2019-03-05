
function [OUT,intermediates]=foot_placement_model_function_step(COM,Rfoot,Lfoot,events,fsopto,pred_samples,order,removeorigin,centerdata)
% function to calculate foot placement model and % explained variance.
% returns output structure with variables and suggested figure titles as
% well as y-labels



% include; -intermediates, n# of steps excluded, steptimes
%% actual calculations. Start with setting things up.
% check the incoming events.
[events,flag] = order_events(events);
if flag~=0
    error('there is something wrong with your events, returning')
end
lhs     = events.lhs;
rhs     = events.rhs;
lto     = events.lto;
rto     = events.rto;

% first find some strides that need to be excluded (any that are too long)
st_L = diff(lhs);
st_R = diff(rhs);


OUT.stride_time.data        = nanmean(st_L./fsopto);
OUT.stride_time_var.data    = nanstd(st_L./fsopto);

OUT.stride_time.titel       ='Stride time';
OUT.stride_time.ylabel      ='Stride time [s]';

OUT.stride_time_var.titel   ='Stride time variability';
OUT.stride_time_var.ylabel  ='Stride time variability [s]';

% calculate velocities
COM_vel = calc_derivative(COM,fsopto);
COM_acc = calc_derivative(COM_vel,fsopto);

% time normalize; left variables
% left swing is lto(1)-lhs(2)
% left stance is rto(2)-rhs(2)
COM_L       = Normalizetimebase_step(COM,lto(1:end-1),lhs(2:end));
COM_L_vel   = Normalizetimebase_step(COM_vel,lto(1:end-1),lhs(2:end));
COM_L_acc   = Normalizetimebase_step(COM_acc,lto(1:end-1),lhs(2:end));
foot_L      = Normalizetimebase_step(Lfoot,rto(2:end),rhs(2:end));
origin_L    = Normalizetimebase_step(Rfoot,lto(1:end-1),lhs(2:end));

% do the above two things also for right
% right swing is rto(1)-rhs(1)
% right stance is lto(1)-lhs(2)
COM_R       = Normalizetimebase_step(COM,rto(1:end-1),rhs(1:end-1));
COM_R_vel   = Normalizetimebase_step(COM_vel,rto(1:end-1),rhs(1:end-1));
COM_R_acc   = Normalizetimebase_step(COM_acc,rto(1:end-1),rhs(1:end-1));
foot_R      = Normalizetimebase_step(Rfoot,lto(1:end-1),lhs(2:end));
origin_R    = Normalizetimebase_step(Lfoot,rto(1:end-1),rhs(1:end-1));

% remove last two colums, and any strides that were too long/short.
COM_L(:,end-1:end)       = nan;
COM_L_vel(:,end-1:end)   = nan;
COM_L_acc(:,end-1:end)   = nan;
foot_L(:,end-1:end)      = nan;
origin_L(:,end-1:end)    = nan;
% remove last two colums
COM_R(:,end-1:end)       = nan;
COM_R_vel(:,end-1:end)   = nan;
COM_R_acc(:,end-1:end)   = nan;
foot_R(:,end-1:end)      = nan;
origin_R(:,end-1:end)    = nan;
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
foot_L       = foot_L(25,:)';
origin_L     = origin_L(25,:)';
foot_R       = foot_R(25,:)';
origin_R     = origin_R(25,:)';

if removeorigin
    foot_L       = foot_L -origin_L; % substract origin, which is the other foot
    foot_R       = foot_R -origin_R ;
end
OUT.SW.data         = nanmean(abs([foot_L; foot_R]));
OUT.SW_var.data     = nanstd(abs([foot_L ;foot_R]));

OUT.SW.ylabel       ='Stepwidth [m]';
OUT.SW_var.ylabel   ='Stepwidth variability[m]';
OUT.SW.titel        ='Stepwidth';
OUT.SW_var.titel    ='Stepwidth variability';

if centerdata
    foot_L       = foot_L -nanmean(foot_L); % substract mean
    foot_R       = foot_R -nanmean(foot_R);
end

%%
L_jac       = zeros(length(pred_samples),order);
R_jac       = zeros(length(pred_samples),order);
warning off
for i_pred_sample=pred_samples
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
    end
    if centerdata
        pred_Lstance        = pred_Lstance-repmat(nanmean(pred_Lstance),size(pred_Lstance,1),1);
        pred_Rstance        = pred_Rstance-repmat(nanmean(pred_Rstance),size(pred_Rstance,1),1);
    end
    foot_L_sample=foot_L;
    foot_R_sample=foot_R;
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
    OUT.Right_N.data(i_pred_sample)      = size(pred_Lstance,1);
    
    stats = regstats(foot_R_sample,pred_Rstance(:,1:order),'linear');
    OUT.Right_pct.data(i_pred_sample)            = stats.rsquare;
    OUT.Right_pct.pval(i_pred_sample)            = stats.fstat.pval;
    OUT.Right_coeff1.data(i_pred_sample)         = stats.beta(2);
    OUT.Right_coeff2.data(i_pred_sample)         = stats.beta(3);
    OUT.Right_coeff1.pval(i_pred_sample)         = stats.tstat.pval(2);
    OUT.Right_coeff2.pval(i_pred_sample)         = stats.tstat.pval(3);
    intermediates.error_right(i_pred_sample,ind_R)   = stats.r;
    
    
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
    
    %% for left
    tmp                                 = [foot_L_sample pred_Lstance];
    ind_L                               = 1:size(tmp,1);
    foot_L_sample(isnan(sum(tmp,2)))    = [];
    pred_Lstance(isnan(sum(tmp,2)),:)   = [];
    ind_L(isnan(sum(tmp,2)))            = [];
    OUT.Left_N.data(i_pred_sample)      = size(pred_Lstance,1);
    
    stats = regstats(foot_L_sample,pred_Lstance(:,1:order),'linear');
    OUT.Left_pct.data(i_pred_sample)            = stats.rsquare;
    OUT.Left_pct.pval(i_pred_sample)            = stats.fstat.pval;
    OUT.Left_coeff1.data(i_pred_sample)         = stats.beta(2);
    OUT.Left_coeff2.data(i_pred_sample)         = stats.beta(3);
    OUT.Left_coeff1.pval(i_pred_sample)         = stats.tstat.pval(2);
    OUT.Left_coeff2.pval(i_pred_sample)         = stats.tstat.pval(3);
    intermediates.error_left(i_pred_sample,ind_L)   = stats.r;
    
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
    %% combined
    foot_combined                       = [foot_R_sample ;foot_L_sample];
    pred_combined                       = [pred_Rstance(:,1:order); pred_Lstance(:,1:order)];
    tmp                                 = [foot_combined pred_combined];
    ind_combined                        = 1:size(tmp,1);
    foot_combined(isnan(sum(tmp,2)))    = [];
    pred_combined(isnan(sum(tmp,2)),:)  = [];
    ind_combined(isnan(sum(tmp,2)))     = [];
    OUT.Combined_N.data(i_pred_sample)      = size(pred_Lstance,1);
    
    
    stats = regstats(foot_combined,pred_combined(:,1:order),'linear');
    OUT.Combined_pct.data(i_pred_sample)            = stats.rsquare;
    OUT.Combined_pct.pval(i_pred_sample)            = stats.fstat.pval;
    OUT.Combined_coeff1.data(i_pred_sample)         = stats.beta(2);
    OUT.Combined_coeff2.data(i_pred_sample)         = stats.beta(3);
    OUT.Combined_coeff1.pval(i_pred_sample)         = stats.tstat.pval(2);
    OUT.Combined_coeff2.pval(i_pred_sample)         = stats.tstat.pval(3);
    intermediates.error_combined(i_pred_sample,ind_combined)   = stats.r;
    
    
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
end
warning off
intermediates.error_right(intermediates.error_right==0)=nan;
intermediates.error_left(intermediates.error_left==0)=nan;
intermediates.error_combined(intermediates.error_combined==0)=nan;

