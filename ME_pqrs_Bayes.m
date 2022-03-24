function [ab_post, p_lab, LME, k] = ME_pqrs_Bayes(y, x, m, ab_prior)
% _
% Parameter estimation for pqrs model using posterior distributions
% FORMAT [ab_post, p_lab, LME, k] = ME_pqrs_Bayes(y, x, m, ab_prior)
% 
%     y        - an n x 1 vector of subsequent memory confidence ratings (1-5)
%     x        - an n x 1 vector of 1s ("old stimuli") and 2s ("new stimuli")
%     m        - a string specifying the computational model (see below)
%     ab_prior - a 1 x 2 vector with prior distribution parameters
%     
%     ab_post  - a k x 2 matrix with posterior distribution parameters
%     p_lab    - a k x 1 cell array indicating names of model parameters
%     LME      - a scalar, the log model evidence of the model
%     k        - the number of free model parameters
% 
% [ab_post, p_lab, LME, k] = ME_pqrs_Bayes(y, x, m, ab_prior) estimates the
% pqrs model m for behavioral responses y with item categories x and
% returns posterior distribution parameters ab_post as well as the log
% model evidence LME.
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
% First edit: 23/02/2021, 15:46
%  Last edit: 24/03/2022, 13:13


% Set input parameters to default values
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(m)
    m = 'pqqrrss';
end;
if nargin < 4 || isempty(ab_prior)
    ab_prior = [1,1];
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

% Calculate posterior distribution parameters
%-------------------------------------------------------------------------%
ab_post = [];
p_lab   = [];
% p - probability of a non-neutral response
if numel(strfind(m,'p')) == 1
    ab_post = [ab_post; ab_prior+[(o1245+n1245), (o3+n3)]];
    p_lab   = [p_lab; {'p'}];
elseif numel(strfind(m,'p')) == 2
    ab_post = [ab_post; ab_prior+[(o1245), (o3)]; ab_prior+[(n1245), (n3)]];
    p_lab   = [p_lab; {'p_o'; 'p_n'}];
else
    warning('p is misspecified (is allowed one or two times)!');
end;
% q - probability of an affirmative, given non-neutral response
if numel(strfind(m,'q')) == 1
    ab_post = [ab_post; ab_prior+[(o45+n45), (o12+n12)]];
    p_lab   = [p_lab; {'q'}];
elseif numel(strfind(m,'q')) == 2
    if numel(strfind(m,'-')) == 0
        ab_post = [ab_post; ab_prior+[(o45), (o12)]; ab_prior+[(n45), (n12)]];
        p_lab   = [p_lab; {'q_o'; 'q_n'}];
    elseif numel(strfind(m,'-')) == 1
        ab_post = [ab_post; ab_prior+[(o45+n12), (o12+n45)]];
        p_lab   = [p_lab; {'q_o=1-q_n'}];
    else
        warning('- is misspecified (is allowed zero or once)!');
    end;
else
    warning('q is misspecified (is allowed one or two times)!');
end;
% r/s - probability of a confident, given affirmative/non-affirmative response
if numel(strfind(m,'r')) == 1
    if numel(strfind(m,'s')) == 0
        ab_post = [ab_post; ab_prior+[(o15+n15), (o24+n24)]];
        p_lab   = [p_lab; {'r'}];
    elseif numel(strfind(m,'s')) == 1
        ab_post = [ab_post; ab_prior+[(o5+n5), (o4+n4)]; ab_prior+[(o1+n1), (o2+n2)]];
        p_lab   = [p_lab; {'r'; 's'}];
    else
        warning('s is misspecified relative to r (the count of s has to match the count of r)!');
    end;
elseif numel(strfind(m,'r')) == 2
    if numel(strfind(m,'s')) == 0
        ab_post = [ab_post; ab_prior+[(o15), (o24)]; ab_prior+[(n15), (n24)]];
        p_lab   = [p_lab; {'r_o'; 'r_n'}];
    elseif numel(strfind(m,'s')) == 2
        if numel(strfind(m,'=')) == 0
            ab_post = [ab_post; ab_prior+[(o5), (o4)]; ab_prior+[(n5), (n4)]; ab_prior+[(o1), (o2)]; ab_prior+[(n1), (n2)]];
            p_lab   = [p_lab; {'r_o'; 'r_n'; 's_o'; 's_n'}];
        elseif numel(strfind(m,'=')) == 1
            if numel(strfind(m,'=_')) == 1
                ab_post = [ab_post; ab_prior+[(o5+n1), (o4+n2)]; ab_prior+[(n5), (n4)]; ab_prior+[(o1), (o2)]];
                p_lab   = [p_lab; {'r_o=s_n'; 'r_n'; 's_o'}];
            elseif numel(strfind(m,'_=')) == 1
                ab_post = [ab_post; ab_prior+[(o5), (o4)]; ab_prior+[(n5+o1), (n4+o2)]; ab_prior+[(n1), (n2)]];
                p_lab   = [p_lab; {'r_o'; 'r_n=s_o'; 's_n'}];
            else
                warning('= is misspecified relative to _ (= must either come before or after _)!');
            end;
        elseif numel(strfind(m,'=')) == 2
            ab_post = [ab_post; ab_prior+[(o5+n1), (o4+n2)]; ab_prior+[(n5+o1), (n4+o2)];];
            p_lab   = [p_lab; {'r_o=s_n'; 'r_n=s_o'}];
        else
            warning('= is misspecified (is allowed zero, one or two times)!');
        end;
    else
        warning('s is misspecified relative to r (the count of s has to match the count of r)!');
    end;
else
    warning('r is misspecified (is allowed one or two times)!');
end;

% Calculate log model evidence
%-------------------------------------------------------------------------%
k   = size(ab_post,1);
LME = 0;
for j = 1:k
    LME = LME + betaln(ab_post(j,1),ab_post(j,2)) - betaln(ab_prior(1),ab_prior(2));
end;