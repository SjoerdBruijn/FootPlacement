clear all;
addpath(genpath('plots'))
load 'testdata'

pred_samples = 1:51;
order        = 2;
removeorigin = 1;
centerdata   = 1;

[OUT,intermediates]=foot_placement_model_function_step(CoM_ML,Rfoot,Lfoot,events,fs_opto,pred_samples,order,removeorigin,centerdata)

figure;
plot((1:51)*2-2,OUT.Combined_pct.data*100)
ylabel(OUT.Combined_pct.ylabel)
title(OUT.Combined_pct.titel)
xlabel('step percentage (%)')

%% create some nice animated plots
FPmodelTimeseriesPlot(CoM_ML,Lfoot,Rfoot,events,fs_opto)


FPmodelAnimatedPlot(CoM_ML,Rfoot,Lfoot,events,fs_opto,'test')