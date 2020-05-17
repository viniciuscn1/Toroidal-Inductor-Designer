function fit=tfid_fit(parameters,data)
% tffit calculates fitness of transfer function fit
%
% t=tf_fit(parameters,data)
%
% Inputs:
% parameters = paramameter vector
%              (first gains, then time constants)
% data       = data.s - complex frequency vector
%              data.t - transfer function values
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
% partial fraction transfer function fitness

o=(length(parameters))/2;

a   = parameters(1:o);
tau = parameters(o+1:2*o);

tpred = pftf(data.s,a,tau);

terror = data.t-tpred;

error = norm(terror./abs(data.t))/length(tpred);

fit = 1.0/(1.0e-12+error);
