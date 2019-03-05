function  [Cycle,TimeGain]=normalizeTimeBase_step(signal,to,hs)
% NormalizeTimeBase: calcultes a 0-100% timebase, ensemble-averages cyclic signals
% [Cycle,TimeGain]=NormalizeTimeBase(signal,trigind)
%
% Input : Signal: any one-dimensional array
%         trigind : array of indices, default: [1 length(Signal)]
%                   should increase monotoneusly
% Process: calculates new points based on a 0-100% time base
%          by spline interpolation for each time interval
% Output:  if length(trigind)=2: Cycle [101 1]
%          if length(trigind)>2: Cycle [101 Ncycles+2]
%             Ncyles=length(trigind)-1,
%             Ncycles+1: mean signal per point, i.e. ensemble averaged
%             Ncycles+2: stand.dev ensemble averaged points
%          TimeGain: (average) amplification/reduction of time-axis (i.e. 100/(samples/cycle))
%
% WARNING user should be aware of information loss in case of excessive downsampling

% AUTHOR(S) AND VERSION-HISTORY
% Ver 1.2 April 2003 (Jaap Harlaar VUmc Amsterdam) adapted from some version
% Ver 1.3 December 2018 (Moira van Leeuwen and Mohammad Mahaki
% Amsterdam) adapted for step 0-100%

% Making sure index first toe off precedes heelstrike
if to(1)>hs(1)
   to = to(1:length(to)-1);
   hs = hs(2:length(hs));
end

% When to and hs are not of equal length
if length(to) > length(hs)
    to = to(1:length(hs));
end

if length(hs) > length(to)
    hs = hs(1:length(to));
end

% if nargin < 2, trigind=[1 length(signal)]; end
% if nargin < 1, return, end

% FFE  a check for validity of indices
nansignal(isnan(signal))=0;
nansignal(~isnan(signal))=1;


Cycle=[1:51]'*nan;
Cyclenan=[1:51]'*nan;
CycleLength=-51;
N=length(to);
if N>1
    for i= 1:N
        x=[to(i):hs(i)]-to(i);
        CycleLength(i)=length(x);
        x=x*50/(hs(i)-to(i));
        x=x+1;
        try
            Cycle(:,i)=interp1(x',signal(to(i):hs(i))',[1:51]','pchip');
            Cyclenan(:,i)=interp1(x',nansignal(to(i):hs(i))',[1:51]','pchip');
            
        catch
            Cyclenan(:,i)=nan;
            Cycle(:,i)=nan;
        end
    end
    Cyclenan(Cyclenan<1)=nan;
    Cycle=Cycle.*Cyclenan;
    
    tmp=nanmean(Cycle(:,1:N)');
    Cycle(:,N+1)=tmp';
    tmp=nanstd(Cycle(:,1:N)');
    Cycle(:,N+2)=tmp';
    TimeGain=51/mean(CycleLength);
elseif N==1,
    x=[to(1):hs(1)]-to(1);
    CycleLength=length(x);
    x=x.*50/(hs(1)-to(1))+1;
    try
        Cycle(:,1)=interp1(x',signal(to(1):hs(2))',[1:51]','cubic');
    catch
        Cycle(1:51,1)=nan;
    end
    TimeGain=51/CycleLength;
end
return
% ============================================
 %END % ### NormalizeTimeBase