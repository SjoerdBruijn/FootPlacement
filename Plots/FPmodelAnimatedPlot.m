
function FPmodelAnimatedPlot(com_timseries,Rfoot_timseries,Lfoot_timseries,events,fs,filename)
if nargin==5
    filename=[];
end
%assumes that signals are ordered ML AP (vt)
if size(com_timseries,2)==1
    warning('Assuming only ML timeseries was inputted, setting AP to random position')
    com_timseries=[com_timseries,zeros(size(com_timseries))+0.5];
    Rfoot_timseries=[Rfoot_timseries,zeros(size(Rfoot_timseries))+1];
    Lfoot_timseries=[Lfoot_timseries,zeros(size(Lfoot_timseries))];
end
%% add some paths
addpath(genpath('gif'))
addpath('..')
%% some settings for foot placement model
pred_samples    = 1:51;
order           = 2;
removeorigin    = 1;
centerdata      = 1;
%%
com     = com_timseries(events.rhs,[1 2]);
rfoot   = Rfoot_timseries(events.rhs,:);
lfoot   = Lfoot_timseries(events.rhs,:);

%% if we were to make a function, should start here
% rereference everything to mean left foot position, just to get nicer position in space,
com     = com-mean(lfoot);
rfoot   = rfoot-mean(lfoot);
lfoot   = lfoot-mean(lfoot);

data=[com;lfoot;rfoot];
% initial plot
figure
subplot(1,2,1)
xlims=[min(data(:,1))-abs(0.1*range(data(:,1))) max(data(:,1))+0.1*range(data(:,1))];
ylims=[min(data(:,2))-abs(0.1*range(data(:,2))) max(data(:,2))+0.1*range(data(:,2))];
h1=plot(lfoot(:,1),lfoot(:,2),'go');hold on
h2=plot(rfoot(:,1),rfoot(:,2),'ro');hold on
h3=plot(com(:,1),com(:,2),'ko');hold on
axis equal
set(gca,'box','off','linewidth',2,'xlim',xlims,'Ylim',ylims,'fontsize',12); %SMB; add Ytick and Xtick so that makes sens? 
ylabel('Y-position (m)','fontsize',14)
xlabel('X-position (m)','fontsize',14)
set(gcf,'color',[1 1 1])
pause
if ~isempty(filename)
    gif([filename,'part1.gif'])
    gif([filename,'part2.gif'])
end
% now, lets visualise re-referencing to the left foot.
for i=1:10
    h1.XData=lfoot(:,1)-i*0.1*lfoot(:,1);
    h1.YData=lfoot(:,2)-i*0.1*lfoot(:,2);
    h2.XData=rfoot(:,1)-i*0.1*lfoot(:,1);
    h2.YData=rfoot(:,2)-i*0.1*lfoot(:,2);
    h3.XData=com(:,1)-i*0.1*lfoot(:,1);
    h3.YData=com(:,2)-i*0.1*lfoot(:,2);
    drawnow
    if ~isempty(filename)
        gif
    end
end

pause
if ~isempty(filename)
    gif([filename,'part3.gif'])
end
% zoom in on the plot
data=[h2.XData' h2.YData'; h3.XData' h3.YData'];
des_xlim=[min(data(:,1))-abs(0.1*range(data(:,1))) max(data(:,1))+0.1*range(data(:,1))];
des_ylim=[min(data(:,2))-abs(0.1*range(data(:,2))) max(data(:,2))+0.1*range(data(:,2))];
for i=1:10
    set(gca,'Xlim',[xlims-0.1*i*xlims+0.1*i*des_xlim])
    set(gca,'Ylim',[ylims-0.1*i*ylims+0.1*i*des_ylim])
    drawnow
    if ~isempty(filename)
        gif
    end
end
% create new variables that are also re-reference (could also get
% the data from the figure, would be nicer
% rfoot2  = rfoot-lfoot;
% com2    = com-lfoot;
% lfoot2  = lfoot-lfoot;
% pause
% if ~isempty(filename)
%     gif([filename,'part4.gif'])
% end
% % 'squash' the figure, to show we are only looking at ML prediction
% here
% for i=1:10
%     h2.YData=rfoot2(:,2)-i*0.1*rfoot2(:,2)+i*0.1*mean(rfoot2(:,2));
%     h3.YData=com2(:,2)-i*0.1*com2(:,2)+i*0.1*mean(com2(:,2));
%     drawnow
%     if ~isempty(filename)
%         gif
%     end
% end

pause
if ~isempty(filename)
    gif([filename,'part5.gif'],'DelayTime',1/50)
end
[OUT,intermediates]=foot_placement_model_function_step(com_timseries,Rfoot_timseries,Lfoot_timseries,events,fs,pred_samples,order,removeorigin,centerdata);
clear data; 
data(:,1)=rfoot(1:end-1,1);
data(:,2)=rfoot(1:end-1,1)+intermediates.error_right(end,:)';
des_xlim=[min(data(:,1))-abs(0.1*range(data(:,1))) max(data(:,1))+0.1*range(data(:,1))];
des_ylim=[min(data(:,2))-abs(0.1*range(data(:,2))) max(data(:,2))+0.1*range(data(:,2))];
% show connection between CoM and foot placement
for i=1:length(com)-1
    subplot(1,2,1)
    plot([h2.XData(i),h3.XData(i)]',[h2.YData(i),h3.YData(i)]','b','linewidth',0.5);hold on
    subplot(1,2,2)
    plot(rfoot(i,1),rfoot(i,1)+intermediates.error_right(end,i),'b.'); hold on
    if i==1
        axis equal
        set(gca,'box','off','linewidth',2,'xlim',des_xlim,'Xtick',[ ],'Ylim',des_ylim,'Ytick',[],...
            'XColor','R','fontsize',14)
        xlabel('Actual foot placement','fontsize',14)
        ylabel('Predicted foot placement','fontsize',14)
    end
    drawnow
    if ~isempty(filename)
        gif
    end
end

Ps = polyfit(rfoot(1:end-1,1),rfoot(1:end-1,1)'+intermediates.error_right(end,:),1);
Ys = polyval(Ps,min(rfoot(:,1)):0.01:max(rfoot(:,1))+0.01);
plot(min(rfoot(:,1)):0.01:max(rfoot(:,1))+0.01,Ys,'b','linewidth',3);
text(des_xlim(1)+0.05,des_ylim(2),sprintf('R^2=%0.2f',OUT.Right_pct.data(end)),'Fontsize',14)
if ~isempty(filename)
    gif
end


