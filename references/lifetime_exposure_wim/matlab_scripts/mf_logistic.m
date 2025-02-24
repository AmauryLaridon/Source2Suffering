% Logistic: calculates the logistic function of the input
% by Will Dwinnell
%
% Last modified: Sep-02-2006

function Output = mf_logistic(Input)

Output = 1 ./ (1 + exp(-Input));

end