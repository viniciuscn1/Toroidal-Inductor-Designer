function [GAP,GAS,Pk]=gainit(objhandle,D,GAP,GAS,iP)
% GAINIT  initialize the genetic algorithm
% [GAP,GAS,Pk] = gainit(@fitfun,D,GAP,GAS,iP)
%
% Inputs:
% fitfun   = name of *.m function that evaluates the fitness.  This
%            function can have up to three vector arguments, one for
%            integer variables, one for real variables which vary over
%            a linear range, and one for exponential variables that
%            vary over a exponential range.  The output should be a simple
%            scalar with the evaluated objective function
% D        = data needed by fitness function (this could be an empty matrix)
% GAP      = Genetic Algorithm Parameters
%            (see gap_default for definitions)
% GAS      = Genetic Algorithm Statistics
% iP       = Initial population (could be empty)
%
% Outputs:
% GAP      = Genetic Algorithm Parameters
% GAS      = Genetic Algorithm Statistics
% Pk       = Population
%
% Written by:
% S.D. Sudhoff
% Purdue University
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize the random number generators
rand('state',sum(100*clock));
randn('state',sum(100*clock));

   
% initialize objective weighting vector
GAP.owv = objwght(GAP);

% set up the population
if isempty(iP)
   [Pk,GAS]=gasetup(GAP.fp_ipop,GAP,objhandle,D);
      else
   Pk=iP;
   if isempty(D)
      [Pk.mfit,Pk.eval,GAS.ne(GAS.cg)]=evaluate(Pk,GAP,GAS.ne(GAS.cg));
   else
      [Pk.mfit,Pk.eval,GAS.ne(GAS.cg)]=evaluate(Pk,GAP,GAS.ne(GAS.cg),D);
   end
end
   
% reduce population size if appropriate
if (GAP.fp_npop < Pk.size)
   Pk2=downsize(Pk,GAP.fp_npop);
   clear Pk;
   Pk=Pk2;
end
    
% report on initial evaluation
GAS=updatestat(GAP,GAS,Pk);
reportplot(GAP,GAS,Pk);
GAS.cg=GAS.cg+1;
       
end

%  Copyright 2013 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of GOSET
% 
%  GOSET is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  GOSET is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with GOSET.  If not, see <http://www.gnu.org/licenses/>.