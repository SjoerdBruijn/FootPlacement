function FPmodelTimeseriesPlot(COM,Lfoot,Rfoot,events,fs)



COM=COM(events.lhs(1):end);
Lfoot=Lfoot(events.lhs(1):end);
Rfoot=Rfoot(events.lhs(1):end);

events.lto=events.lto-events.lhs(1)+1;
events.rto=events.rto-events.lhs(1)+1;
events.rhs=events.rhs-events.lhs(1)+1;
events.lhs=events.lhs-events.lhs(1)+1;
%%
addpath('..')
%% some settings for foot placement model
pred_samples    = 1:51;
order           = 2;
removeorigin    = 1;
centerdata      = 1;



ms_l =round((events.lhs+events.lto)/2);
ms_r =round((events.rhs(1:end-1)+events.rto(2:end))/2);
%select only first 10 seconds



y_l = Lfoot(ms_l);
y_r = Rfoot(ms_r);
com_y_l =COM(ms_l);
com_y_r =COM(ms_r);
plot((1:length(COM))/fs,COM,'linewidth',2,'Color',[0 0 1]);hold on % com plot
plot(ms_l/fs,com_y_l,'o','Color',[0 0.7 0],'LineWidth',2) % CoM at left midstance
plot(ms_r/fs,com_y_r,'o','Color',[1 0 0],'LineWidth',2) % right
plot([events.lhs]'./fs,[Lfoot(events.lhs)]','o','Color',[0 0.7 0],'MarkerFaceColor', [0 0.7 0],'linewidth',2); % left foot
plot(events.rhs(1:end-1)'./fs,Rfoot(events.rhs(1:end-1))','ro','MarkerFaceColor', [1 0 0],'linewidth',2) % right foot



[~,intermediates]=foot_placement_model_function_step(COM,Rfoot,Lfoot,events,2,pred_samples,order,removeorigin,centerdata);
y_l_err=Lfoot(events.lhs(2:end))+intermediates.error_left(25,:)';
y_r_err=Rfoot(events.rhs(2:end))+intermediates.error_right(25,:)';


plot(events.lhs(2:end)/fs ,y_l_err,'g*')%,'Color',[0 0.7 0]) % left foot predictions
plot(events.rhs(1:end-1)/fs ,y_r_err,'r*') % right foot predictions
axis('manual')
set(gca,'xlim',[0.5 8.5], ...
    'Xtick',[0.5:1:8.8],...
    'Xticklabel',[0:1:8],...
    'PlotBoxAspectRatio',[5,1.5,1], ...
    'box','off', ...
    'LineWidth',2 ,...
    'fontsize',14, ...
    'YTick',[])

arrow3([ms_l(1:end-1)/fs com_y_l(1:end-1)],[ events.rhs(1:end-1)/fs y_r_err ])
arrow3([ms_r(1:end)/fs com_y_r(1:end)],[ events.lhs(2:end)/fs y_l_err ])

xlabel('Time (s)','fontsize',14)
ylabel('Position (m)','fontsize',14)


