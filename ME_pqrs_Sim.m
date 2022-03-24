function [y, x, m] = ME_pqrs_Sim(p, p_lab, n, o)
% _
% Simulate behavioral responses using pqrs model
% FORMAT [y, x, m] = ME_pqrs_Sim(p, p_lab, n, o)
% 
%     p     - a k x 1 vector of model parameter values
%     p_lab - a k x 1 cell array indicating names of model parameters
%     n     - a integer, the number of trials to simulate (e.g. 90)
%     o     - a scalar, the fraction of old items (e.g. 2/3)
% 
%     y     - an n x 1 vector of subsequent memory confidence ratings (1-5)
%     x     - an n x 1 vector of 1s ("old stimuli") and 2s ("new stimuli")
%     m     - a string specifying the computational model (see below)
% 
% [y, x, m] = ME_pqrs_Sim(p, p_lab, n, o) simulates item statuses x and
% behavioral responses y for n trials using parameter values p and
% parameter labels p_lab implying model m, assuming old item frequency o.
% 
% NOTE: The product of o and n must be a whole number!
% 
% Author: Joram Soch, DZNE Göttingen
% E-Mail: Joram.Soch@DZNE.de
% 
% First edit: 19/11/2021, 11:29
%  Last edit: 24/03/2022, 13:11


% Create parameter vector
%-------------------------------------------------------------------------%
pi    = p;                      % parameters from input
li    = p_lab;                  % labels from input
p_lab = {'p_o', 'p_n', 'q_o', 'q_n', 'r_o', 'r_n', 's_o', 's_n'}';
p     = NaN(size(p_lab));       % parsed parameter values
m     = '';                     % parsed model name

% Parse parameter values
%-------------------------------------------------------------------------%
for j = 1:numel(li)
    if     strcmp(li{j},'p')                % p
        p([1,2]) = pi(j);
        m = strcat(m,'p');
    elseif strcmp(li{j},'p_o')
        p(1) = pi(j);
        m = strcat(m,'p');
    elseif strcmp(li{j},'p_n')
        p(2) = pi(j);
        m = strcat(m,'p');
    elseif strcmp(li{j},'q')                % q
        p([3,4]) = pi(j);
        m = strcat(m,'q');
    elseif strcmp(li{j},'q_o')
        p(3) = pi(j);
        m = strcat(m,'q');
    elseif strcmp(li{j},'q_n')
        p(4) = pi(j);
        m = strcat(m,'q');
    elseif strcmp(li{j},'q_o=1-q_n')        % qq-
        p(3) = pi(j);
        p(4) = 1-pi(j);
        m = strcat(m,'qq-');
    elseif strcmp(li{j},'r')                % r
        p([5,6]) = pi(j);
        m = strcat(m,'r');
    elseif strcmp(li{j},'r_o')
        p(5) = pi(j);
        m = strcat(m,'r');
    elseif strcmp(li{j},'r_n')
        p(6) = pi(j);
        m = strcat(m,'r');
    elseif strcmp(li{j},'s')                % s
        p([7,8]) = pi(j);
        m = strcat(m,'s');
    elseif strcmp(li{j},'s_o')
        p(7) = pi(j);
        m = strcat(m,'s');
    elseif strcmp(li{j},'s_n')
        p(8) = pi(j);
        m = strcat(m,'s');
    elseif strcmp(li{j},'r_o=s_n')          % rrss==
        p([5,8]) = pi(j);
        m = strcat(m,'=');
    elseif strcmp(li{j},'r_n=s_o')
        p([6,7]) = pi(j);
        m = strcat(m,'=');
    else
        warning('Parameter "%s" is misspecified (label doesn''t exist)!', li{j});
    end;
end;
if isnan(p(7)), p(7) = p(5); end;           % no s, only r
if isnan(p(8)), p(8) = p(6); end;

% Simulate item statuses
%-------------------------------------------------------------------------%
X = [[1*ones(round(o*n),1); 2*ones(n-round(o*n),1)], rand(n,1)];
X = sortrows(X,2);
x = X(:,1);                     % item status: 1 = old, 2 = new
clear X

% Simulate behavioral responses
%-------------------------------------------------------------------------%
R = rand(n,3);                  % decided, affirmative, confident
y = zeros(n,1);                 % behavioral response: 1-5
y(x==1 & R(:,1)< p(1) & R(:,2)>=p(3) & R(:,3)< p(7)) = 1;
y(x==1 & R(:,1)< p(1) & R(:,2)>=p(3) & R(:,3)>=p(7)) = 2;
y(x==1 & R(:,1)>=p(1))                               = 3;
y(x==1 & R(:,1)< p(1) & R(:,2)< p(3) & R(:,3)>=p(5)) = 4;
y(x==1 & R(:,1)< p(1) & R(:,2)< p(3) & R(:,3)< p(5)) = 5;
y(x==2 & R(:,1)< p(2) & R(:,2)>=p(4) & R(:,3)< p(8)) = 1;
y(x==2 & R(:,1)< p(2) & R(:,2)>=p(4) & R(:,3)>=p(8)) = 2;
y(x==2 & R(:,1)>=p(2))                               = 3;
y(x==2 & R(:,1)< p(2) & R(:,2)< p(4) & R(:,3)>=p(6)) = 4;
y(x==2 & R(:,1)< p(2) & R(:,2)< p(4) & R(:,3)< p(6)) = 5;