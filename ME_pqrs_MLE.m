function [p_MLE, p_lab, MLL, k] = ME_pqrs_MLE(y, x, m)
% _
% Parameter estimation for pqrs model using maximum likelihood estimation
% FORMAT [p_MLE, p_lab, MLL, k] = ME_pqrs_MLE(y, x, m)
% 
%     y     - an n x 1 vector of subsequent memory confidence ratings (1-5)
%     x     - an n x 1 vector of 1s ("old stimuli") and 2s ("new stimuli")
%     m     - a string specifying the computational model (see below)
%     
%     p_MLE - a k x 1 vector of maximum likelihood parameter estimates
%     p_lab - a k x 1 cell array indicating names of model parameters
%     MLL   - a scalar, the maximum log-likelihood of the model
%     k     - the number of free model parameters
% 
% [p_MLE, p_lab, MLL, k] = ME_pqrs_MLE(y, x, m) estimates the pqrs model m
% for behavioral responses y with item categories x and returns maximum
% likelihood estimates p_MLE as well as the maxmimum log-likelihood MLL.
% 
% The input variable "m" is a string such as "pqqrs" specifying the model:
% o p is the probability of a non-neutral response
% o q is the probability of an affirmative, given non-neutral response
% o r is the probability of a confident, given affirmative response
% o s is the probability of a confident, given non-affirmative response
% - If a letter is listed twice, the probability is estimated separately
%   for old and new stimuli.
% - If the letter "s" does not occur in the model specification, it is
%   assumed to be equal with "r".
% - If "q" occurs twice, there is the possibility to constrain
%   - q_new = 1-q_old via "-" in m.
% - If "r" and "s" both occur twice, there is the possibility to constrain
%   - r_old = s_new via "=_" in m;
%   - r_new = s_old via "_=" in m;
%   - r_old = s_new and r_new = s_old via "==" in m.
% 
% Author: Joram Soch, DZNE Göttingen
% E-Mail: Joram.Soch@DZNE.de
% 
% First edit: 23/02/2021, 12:54
%  Last edit: 24/03/2022, 13:12


% Set input parameters to default values
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(m)
    m = 'pqqrrss';
end;

% Obtain behavioral summary statistics
%-------------------------------------------------------------------------%
ON = zeros(2,5);
for type = 1:2
    for resp = 1:5
        ON(type,resp) = sum(x==type & y==resp);
    end;
end;
o1 = ON(1,1); o2 = ON(1,2); o3 = ON(1,3); o4 = ON(1,4); o5 = ON(1,5);
n1 = ON(2,1); n2 = ON(2,2); n3 = ON(2,3); n4 = ON(2,4); n5 = ON(2,5);
o12 = o1+o2; o45 = o4+o5; o15 = o1+o5; o24 = o2+o4;
n12 = n1+n2; n45 = n4+n5; n15 = n1+n5; n24 = n2+n4;
o1245 = o12+o45;
n1245 = n12+n45;

% Calculate maximum likelihood estimates
%-------------------------------------------------------------------------%
p_MLE = [];
p_lab = [];
% p - probability of a non-neutral response
if numel(strfind(m,'p')) == 1
    p_MLE = [p_MLE; (o1245+n1245)/(o1245+o3+n1245+n3)];
    p_lab = [p_lab; {'p'}];
elseif numel(strfind(m,'p')) == 2
    p_MLE = [p_MLE; (o1245)/(o1245+o3); (n1245)/(n1245+n3)];
    p_lab = [p_lab; {'p_o'; 'p_n'}];
else
    warning('p is misspecified (is allowed one or two times)!');
end;
% q - probability of an affirmative, given non-neutral response
if numel(strfind(m,'q')) == 1
    p_MLE = [p_MLE; (o45+n45)/(o1245+n1245)];
    p_lab = [p_lab; {'q'}];
elseif numel(strfind(m,'q')) == 2
    if numel(strfind(m,'-')) == 0
        p_MLE = [p_MLE; (o45)/(o1245); (n45)/(n1245)];
        p_lab = [p_lab; {'q_o'; 'q_n'}];
    elseif numel(strfind(m,'-')) == 1
        p_MLE = [p_MLE; (o45+n12)/(o1245+n1245)];
        p_lab = [p_lab; {'q_o=1-q_n'}];
    else
        warning('- is misspecified (is allowed zero or once)!');
    end;
else
    warning('q is misspecified (is allowed one or two times)!');
end;
% r/s - probability of a confident, given affirmative/non-affirmative response
if numel(strfind(m,'r')) == 1
    if numel(strfind(m,'s')) == 0
        p_MLE = [p_MLE; (o15+n15)/(o1245+n1245)];
        p_lab = [p_lab; {'r'}];
    elseif numel(strfind(m,'s')) == 1
        p_MLE = [p_MLE; (o5+n5)/(o45+n45); (o1+n1)/(o12+n12)];
        p_lab = [p_lab; {'r'; 's'}];
    else
        warning('s is misspecified relative to r (the count of s has to match the count of r)!');
    end;
elseif numel(strfind(m,'r')) == 2
    if numel(strfind(m,'s')) == 0
        p_MLE = [p_MLE; (o15)/(o1245); (n15)/(n1245)];
        p_lab = [p_lab; {'r_o'; 'r_n'}];
    elseif numel(strfind(m,'s')) == 2
        if numel(strfind(m,'=')) == 0
            p_MLE = [p_MLE; (o5)/(o45); (n5)/(n45); (o1)/(o12); (n1)/(n12)];
            p_lab = [p_lab; {'r_o'; 'r_n'; 's_o'; 's_n'}];
        elseif numel(strfind(m,'=')) == 1
            if numel(strfind(m,'=_')) == 1
                p_MLE = [p_MLE; (o5+n1)/(o45+n12); (n5)/(n45); (o1)/(o12)];
                p_lab = [p_lab; {'r_o=s_n'; 'r_n'; 's_o'}];
            elseif numel(strfind(m,'_=')) == 1
                p_MLE = [p_MLE; (o5)/(o45); (n5+o1)/(n45+o12); (n1)/(n12)];
                p_lab = [p_lab; {'r_o'; 'r_n=s_o'; 's_n'}];
            else
                warning('= is misspecified relative to _ (= must either come before or after _)!');
            end;
        elseif numel(strfind(m,'=')) == 2
            p_MLE = [p_MLE; (o5+n1)/(o45+n12); (n5+o1)/(n45+o12)];
            p_lab = [p_lab; {'r_o=s_n'; 'r_n=s_o'}];
        else
            warning('= is misspecified (is allowed zero, one or two times)!');
        end;
    else
        warning('s is misspecified relative to r (the count of s has to match the count of r)!');
    end;
else
    warning('r is misspecified (is allowed one or two times)!');
end;

% Calculate maximum log-likelihood
%-------------------------------------------------------------------------%
MLL = 0;
% p - probability of a non-neutral response
if numel(strfind(m,'p')) == 1
    p_est = (o1245+n1245)/(o1245+o3+n1245+n3); p_est = aNaN(p_est);
    MLL   = MLL + (o1245+n1245)*log(p_est) + (o3+n3)*log(1-p_est);
elseif numel(strfind(m,'p')) == 2
    p_est = [(o1245)/(o1245+o3); (n1245)/(n1245+n3)]; p_est = aNaN(p_est);
    MLL   = MLL + (o1245)*log(p_est(1)) + (o3)*log(1-p_est(1)) + ...
                  (n1245)*log(p_est(2)) + (n3)*log(1-p_est(2));
end;
% q - probability of an affirmative, given non-neutral response
if numel(strfind(m,'q')) == 1
    p_est = (o45+n45)/(o1245+n1245); p_est = aNaN(p_est);
    MLL   = MLL + (o45+n45)*log(p_est) + (o12+n12)*log(1-p_est);
elseif numel(strfind(m,'q')) == 2
    if numel(strfind(m,'-')) == 0
        p_est = [(o45)/(o1245); (n45)/(n1245)]; p_est = aNaN(p_est);
        MLL   = MLL + (o45)*log(p_est(1)) + (o12)*log(1-p_est(1)) + ...
                      (n45)*log(p_est(2)) + (n12)*log(1-p_est(2));
    elseif numel(strfind(m,'-')) == 1
        p_est = [(o45+n12)/(o1245+n1245)]; p_est = aNaN(p_est);
        MLL   = MLL + (o45+n12)*log(p_est) + (o12+n45)*log(1-p_est);
    end;
end;
% r/s - probability of a confident, given affirmative/non-affirmative response
if numel(strfind(m,'r')) == 1
    if numel(strfind(m,'s')) == 0
        p_est = (o15+n15)/(o1245+n1245); p_est = aNaN(p_est);
        MLL   = MLL + (o15+n15)*log(p_est) + (o24+n24)*log(1-p_est);
    elseif numel(strfind(m,'s')) == 1
        p_est = [(o5+n5)/(o45+n45); (o1+n1)/(o12+n12)]; p_est = aNaN(p_est);
        MLL   = MLL + (o5+n5)*log(p_est(1)) + (o4+n4)*log(1-p_est(1)) + ...
                      (o1+n1)*log(p_est(2)) + (o2+n2)*log(1-p_est(2));
    end;
elseif numel(strfind(m,'r')) == 2
    if numel(strfind(m,'s')) == 0
        p_est = [(o15)/(o1245); (n15)/(n1245)]; p_est = aNaN(p_est);
        MLL   = MLL + (o15)*log(p_est(1)) + (o24)*log(1-p_est(1)) + ...
                      (n15)*log(p_est(2)) + (n24)*log(1-p_est(2));
    elseif numel(strfind(m,'s')) == 2
        if numel(strfind(m,'=')) == 0
            p_est = [(o5)/(o45); (n5)/(n45); (o1)/(o12); (n1)/(n12)]; p_est = aNaN(p_est);
            MLL   = MLL + (o5)*log(p_est(1)) + (o4)*log(1-p_est(1)) + ...
                          (n5)*log(p_est(2)) + (n4)*log(1-p_est(2)) + ...
                          (o1)*log(p_est(3)) + (o2)*log(1-p_est(3)) + ...
                          (n1)*log(p_est(4)) + (n2)*log(1-p_est(4));
        elseif numel(strfind(m,'=')) == 1
            if numel(strfind(m,'=_')) == 1
                p_est = [(o5+n1)/(o45+n12); (n5)/(n45); (o1)/(o12)]; p_est = aNaN(p_est);
                MLL   = MLL + (o5+n1)*log(p_est(1)) + (o4+n2)*log(1-p_est(1)) + ...
                              (n5)*log(p_est(2))    + (n4)*log(1-p_est(2)) + ...
                              (o1)*log(p_est(3))    + (o2)*log(1-p_est(3));
            elseif numel(strfind(m,'_=')) == 1
                p_est = [(o5)/(o45); (n5+o1)/(n45+o12); (n1)/(n12)]; p_est = aNaN(p_est);
                MLL   = MLL + (o5)*log(p_est(1))    + (o4)*log(1-p_est(1)) + ...
                              (n5+o1)*log(p_est(2)) + (n4+o2)*log(1-p_est(2)) + ...
                              (n1)*log(p_est(3))    + (n2)*log(1-p_est(3));
            end;
        elseif numel(strfind(m,'=')) == 2
            p_est = [(o5+n1)/(o45+n12); (n5+o1)/(n45+o12)]; p_est = aNaN(p_est);
            MLL   = MLL + (o5+n1)*log(p_est(1)) + (o4+n2)*log(1-p_est(1)) + ...
                          (n5+o1)*log(p_est(2)) + (n4+o2)*log(1-p_est(2));
        end;
    end;
end;

% Return number of free parameters
%-------------------------------------------------------------------------%
k = numel(p_MLE);

% Function: avoid not-a-number
%-------------------------------------------------------------------------%
function p = aNaN(p);

% p(isnan(p)) = 0;
p(p==0) = 0+exp(-23);
p(p==1) = 1-exp(-23);

% Explanation: When p is estimated as 0, then 0*log(0) returns NaN, but it
% should return 0 (as in entropy calculation). To enforce this, without
% damaging the calculation of n*log(1), p is set to a very small value.