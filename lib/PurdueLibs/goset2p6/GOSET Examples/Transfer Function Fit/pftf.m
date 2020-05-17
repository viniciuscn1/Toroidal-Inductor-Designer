function t = pf_tf(s,a,tau)
% PF_tf calculates transfer function values for a transfer function
%       with real poles in partial fraction expansion form.
%
% t=pf_tf(s,a,tau)
%
% Inputs:
% s       = vector of complex frequencies
% a       = vector of dc gains
% tau     = vector of time constants
%           (should have same dimension as a)
% Outputs:
% t       = vector of transfer function values
%           (same dimension as s)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

n=length(a);

if n~=length(tau)
   error('alpha and tau vectors do not correspond');
end

t=zeros(size(s));
for i=1:n,
   t=t+a(i)./(tau(i)*s+1);
end