% % By Reza Mahini may 2017
% The input is  "time sample" and our put is "time(ms)"

% ephStart= start of epoch
% ephEnd = end of epoch
% samples = total number of samples in one epoch
% v =time sample (start ts)
% w =time sample (End ts)


function [start_ms,end_ms]=time_conv_ms(samples,v,w,ephStart,ephEnd)

step=samples/(ephEnd-ephStart);

start_ms=v/step+ephStart; % start time (ms)
end_ms=w/step+ephStart; % end time (ms)

end 


