function [RR, RR_ci, FAR, abs_risk_red, rel_risk_red] = mf_riskratio(reference, experiment, alpha, percentages)

% inspired by
% https://www.mathworks.com/matlabcentral/fileexchange/15347-odds

% This function calculates the Risk Ratio and the Odds Ratio (OR) on a 2x2
% input matrix. Both ratios are computed with confidence intervals. If
% confidence interval of OR doesn't encompass the value OR=1, then the
% function computes the Bayesian Credibility Assessment of the test. If the
% test is credible, the function calculates the Association Parameter Phi.
% The association parameter Phi=sqrt(chisquare/N).
% The routine coumputes the Power and, if necessary, the sample sizes needed
% to achieve a power=0.80 using a modified asymptotic normal method with
% continuity correction as described by Hardeo Sahai and Anwer Khurshid in
% Statistics in Medicine, 1996, Vol. 15, Issue 1: 1-21.


% get number of points of the two vectors
nobs_ref = length(reference);
nobs_exp = length(experiment);


% some number needed for CI calculation
Za = -realsqrt(2) * erfcinv(2-alpha);


% loop over percentages
for i=1:length(percentages)
    
    
    % load variable data
    percentage = percentages(i);            % get percentage

    
    % get associated threshold value 
    threshold = prctile(reference, percentage);
    
    
    % get associated probabilities
    prob_ref = 1 - percentage / 100;
    prob_exp = length(find(experiment >= threshold )) ./ nobs_exp;
    p        = [prob_exp; prob_ref];

     
    % construct x (only leeded for CI computation)
    %           X - 2x2 data matrix composed like this: 
    %...........................................extreme..normal
    %                                              ___________
    %experiment                                   |  A  |  B  |
    %                                             |_____|_____|
    %reference                                    |  C  |  D  |
    %                                             |_____|_____|
    x = [nobs_exp*prob_exp   nobs_exp*(1-prob_exp); ...
         nobs_ref*prob_ref   nobs_ref*(1-prob_ref)     ];    
    R = sum(x,2);    % sum of the rows - should be identical to nobs

    
    % % some checks for testing
    % prob_ref_check  = length(find(reference >= threshold )) ./ nobs_ref
    % threshold_check = prctile(experiment, (1-prob_exp)*100) % use this as a check
    % p_check         = x(:,1) ./ R; % probabilities - should be identical to p


    % Risk ratio (RR)
    RR(i)           = prob_exp / prob_ref;  %#ok<*AGROW>
    RR_se(i,1:2)    = realsqrt(sum(1./x(:,1)-1./R));                   % standard error of log(RR)
    RR_ci(i,1:2)    = exp(reallog(RR(i))+([-1 1].*(Za*RR_se(i,1:2)))); % RR confidence interval
    FAR(i)          = 1- ( prob_ref / prob_exp );                      % Fraction of Attributable Risk
    abs_risk_red(i) = abs(prob_ref - prob_exp);                        % absolute risk reduction
    rel_risk_red(i) = abs_risk_red(i) / prob_exp;                      % relative risk reduction
    
    % % some checks for testing
    % fprintf('Risk Ratio: %0.4f<%0.4f<%0.4f\n',rrci(1),rr,rrci(2))
    % if(rrci(1)<=1 && rrci(2)>=1)
    %    disp('Confidence interval encompasses RR=1. None significative association.')
    % end
    % fprintf('Absolute risk reduction: %0.1f%%\n',d*100)
    % fprintf('Relative risk reduction: %0.1f%%\n',rrr*100)
    % disp(' ')

    if reference == 0
        RR           = NaN; 
        RR_ci        = [NaN NaN]; 
        FAR          = NaN; 
        abs_risk_red = NaN;       
        rel_risk_red = NaN;       
    end
    

end
