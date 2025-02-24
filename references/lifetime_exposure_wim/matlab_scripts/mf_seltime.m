
% --------------------------------------------------------------------
% function to select the time frame for the evaluation data
% --------------------------------------------------------------------


function [T_out] = mf_seltime(T_in, time_begin, time_end)


% find first and last index corresponding to time frame 
t1 = find(T_in(:,1) == time_begin(1), 1, 'first');
t2 = find(T_in(:,1) == time_end(1)  , 1, 'last' );


% security: check presence
if isempty(t1); 
    t1 = 1;
end

if isempty(t2); 
    t2 = size(T_in,1);
end


% select frame
T_out = T_in(t1:t2,:);


end
