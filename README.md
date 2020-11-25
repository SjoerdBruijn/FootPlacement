# FootPlacement
Foot placement prediction, based on the work of Wang &amp; Srinivasavan

'foot_placement_model_function_step' is a function to perform a linear regression model. 
It correlates the center of mass kinematic state during swing with subsequent foot placement.

We uploaded example data 'testdata.mat' which can be used to run the example code below.
Here we run the foot placement model and plot the resulting relative explained variance (% explained variance).

Using the foot placement model function requires the Matlab Statistics toolbox and functions from the =VU 3D model= folder.

```matlab
clear all;
load 'testdata'

addpath(genpath('=VU 3D model='))

pred_samples = 1:51;
order        = 2;
removeorigin = 1;
centerdata   = 1;

[OUT,intermediates]=foot_placement_model_function_step(CoM_ML,Rfoot,Lfoot,events,fs_opto,pred_samples,order,removeorigin,centerdata)

figure;
plot((1:51)*2-2,OUT.Combined_pct.data)
ylabel(OUT.Combined_pct.ylabel)
title(OUT.Combined_pct.titel)
xlabel('step percentage (%)')
```
