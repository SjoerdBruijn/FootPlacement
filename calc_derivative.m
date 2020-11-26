function derivative=calc_derivative(signal,samplefrequency,spline_on)
%function derivative=calc_derivative(signal,samplefrequency,spline_on)
% -------------------------------------------------------------------------
% This function calculates the derivative by using diff(signal)*samplefrequency
% 1/2 sample shift is corrected by linear interpolation by devault
% optionally splipne interpolation is performed for possible better result 
% -------------------------------------------------------------------------
% %%Input
% signal            data         
% samplefrequency   sample frequency
% spline_on         optional spline interpolation (spline_on=1)
% %%Output
% derivative derivative of signal 

derivative_shift=diff(signal)*samplefrequency;%1/2 sample shift because diff is the difference between two sampes

% to correct for the 1/2 sample shift we add a nan at the beginning  so
% that sample 1 1/2 is positioned on sample 2. Then we caclulate the
% average beween subsequnt samples that result in another 1/2 a sample 
% shift forward causing sample 2 to end up at sample 2. Finally, the average is
% calculated by taking the shifted derivatives (derivative_shift) and add
% half of the differences between the subsequent samples. Compicated but it
% works ;) Maybe there is a easier solution?

derivative=repmat(nan,size(signal,1),size(signal,2));
derivative(2:end-1,:)= derivative_shift(1:end-1,:)+diff(derivative_shift)/2;


%code below does not work because nanas from other coloms are projected to
%all colloms
%derivative2(:,col_oke)=interp1(t_shift',derivative_shift(:,col_oke),t_original','spline');
% derivative2(:,col_oke)=spline(t_shift,derivative_shift(:,col_oke)',t_original)';

if exist('spline_on','var')
    if spline_on==1
        %spline must be run repeated for each colm because otherwise it copies every nan
        %to the other colms
        t_original=1:size(signal,1);%time of original signal
        t_shift=1.5:size(signal,1)-.5;% time of derivatative (shifter 1/2 sample backwards)
        col_oke=find(sum(~isnan(signal))>2);
        warning off
        for i_col=col_oke
            derivative(:,i_col)=spline(t_shift,derivative_shift(:,i_col)',t_original)';
        end
        warning on
        derivative(isnan(signal))=nan;
        derivative([1 end],:)=nan;
    end
end

% derivative2(isnan(signal))=nan;
% derivative2([1 end],:)=nan;


%figure;plot(derivative);hold on;plot(derivative2,'r.')

% lastwarn
% figure
% plot(t_original, derivative,'b.-');hold on
% plot(t_shift, derivative_shift,'r.-');hold on
%
% lastwarn
% figure
% plot(t_original, derivative2,'b.-');hold on
% plot(t_shift, derivative_shift,'r.-');hold on
%
