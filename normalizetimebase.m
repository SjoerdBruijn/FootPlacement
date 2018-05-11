function  [Cycle,TimeGain]=normalizetimebase(signal,trigind)

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
% new ver, may 2018, SMB


if nargin < 2, trigind=[1 length(signal)]; end
if nargin < 1, return, end

% FFE  a check for validity of indices
nansignal(isnan(signal))=0;
nansignal(~isnan(signal))=1;

N           = length(trigind)-1;
Cycle       = nan(101,N);
Cyclenan    = nan(101,N);
CycleLength = nan;


for i_cycle = 1:N
    x                       = (trigind(i_cycle):trigind(i_cycle+1))-trigind(i_cycle);
    CycleLength(i_cycle)    = length(x);
    x                       = x*100/(trigind(i_cycle+1)-trigind(i_cycle));
    x                       = x+1;
    try
        Cycle(:,i_cycle)    = interp1(x',signal(trigind(i_cycle):trigind(i_cycle+1))',(1:101)','pchip');
        Cyclenan(:,i_cycle) = interp1(x',nansignal(trigind(i_cycle):trigind(i_cycle+1))',(1:101)','pchip');
    catch
        Cyclenan(:,i_cycle) = nan;
        Cycle(:,i_cycle)    = nan;
    end
end
Cyclenan(Cyclenan<1)    = nan;
Cycle                   = Cycle.*Cyclenan;
TimeGain                = 101/mean(CycleLength);

if N>1
    Cycle(:,N+1)    = nanmean(Cycle,2);
    Cycle(:,N+2)    = nanstd(Cycle,[],2);
end

