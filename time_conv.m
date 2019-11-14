% % By Reza Mahini may 2017
% the input is  "time sample" and our put is "time(ms)"
% the form of usage [start_ms,end_ms]=time_conv(v,samples,ephStart,ephEnd)
% ephStart= start of epoch
% ephEnd = end of epoch
% v = requested time samples
% samples = total number of samples in one epoch



function [start_ms,end_ms]=time_conv(v,samples,ephStart,ephEnd)

step=samples/(ephEnd-ephStart);

start_ms=v/step+ephStart; % start time (ms)
% % % end_ms=w/step+ephStart; % end time (ms)

end 


