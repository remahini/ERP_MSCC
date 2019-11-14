% By Reza Mahini may 2017
% the input is time in "ms" the our put is "time sample"

% tempStart= epoch start
% tempEnd= epoch end
% start_ms= known time(ms)
% end_ms=known time (ms)


function [v,w]=time_conv_ts(samples,start_ms,end_ms,ephStart,ephEnd)

step=samples/(ephEnd-ephStart);

v=abs(ephStart-start_ms)*step; % time sample

w=samples-abs(ephEnd-end_ms)*step; % time sample

end