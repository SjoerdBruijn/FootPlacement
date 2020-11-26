# FootPlacement
Foot placement prediction, based on the work of Wang &amp; Srinivasavan (2014), and as used in Mahaki et al (2019), Van Leeuwen et al (2020). 
Feel free to use this code, but please cite one of the above work, or the DOI itself;

'foot_placement_model_function_step' is a function to perform a linear regression model. 
It correlates the center of mass kinematic state during swing with subsequent foot placement.

We uploaded example data 'testdata.mat' which can be used to run the example code below.
Here we run the foot placement model and plot the resulting relative explained variance (% explained variance).

Using the foot placement model function requires the Matlab Statistics toolbox.

```matlab
clear all;
load 'testdata'

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
References
Wang, Yang, and Manoj Srinivasan. "Stepping in the direction of the fall: the next foot placement can be predicted from current upper body state in steady-state walking." Biology letters 10.9 (2014): 20140405.

Mahaki, Mohammadreza, Sjoerd M. Bruijn, and Jaap H. Van Dieën. "The effect of external lateral stabilization on the use of foot placement to control mediolateral stability in walking and running." PeerJ 7 (2019): e7939.

van Leeuwen, Anina Moira, et al. "Active foot placement control ensures stable gait: Effect of constraints on foot placement and ankle moments." bioRxiv (2020).
