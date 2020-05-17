function [x,f] = trimga(GAP,P,D)
% TRIMGA  Uses the Nelder-Mead simplex algorithm to perform a optimization
%         using the best individual from a GA as a starting point.  
%         Gene range constraints are enforced by subtracting infinity
%         from the fitness function when the gene range goes outside of 
%         the prescribed limits.
%
% [x,f] = trimga(GAP,P,D)
% [x,f] = trimga(GAP,P)
%
% Inputs:
% GAP      = Genetic Algorithm Parameters
%            (see gap_default for definitions)
% P        = Population 
% D        = data needed by fitness function 
%            (this could be an empty matrix)
%
% Outputs:
% x        = optimized solution
% f        = optimized fitness function
%
% Written by:
% Ricky Chan and Aaron Cramer for S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

if (GAP.fp_obj==0)
   error('TRIMGA only works with single-objective optimizations');
end

% identify those fields which are fixed
n = length(GAP.gd_type);
intIndex    = find(GAP.gd_type == 1);
nonIntIndex = setdiff(1:n,intIndex);

% get best individual
[oldfitness,bestFitIndex] = max(P.mfit(GAP.fp_obj,:));

% get fixed genes
intValue    = P.gene(intIndex,bestFitIndex);
nonIntValue = P.gene(nonIntIndex,bestFitIndex);

% make anonymous function
if (nargin > 2) 
    df = @(x) -dummy(n,intIndex,nonIntIndex,intValue,x,P,GAP,D);
else
    df = @(x) -dummy(n,intIndex,nonIntIndex,intValue,x,P,GAP);
end

% optimize the free parameters
[nonIntValue,minusf,~,~] = fminsearch(df,nonIntValue);

% assign solution and fitness
x(intIndex)=intValue;
x(nonIntIndex)=nonIntValue;
f=-minusf;

% search for limit violations;
limits=zeros(size(x));
limits(x>GAP.gd_max)=1;
limits(x<GAP.gd_min)=1;

% some reporting
disp(['TRIMGA increased fitness of objective ' ...
      num2str(GAP.fp_obj) ' from ' num2str(oldfitness) ' to ' ...
      num2str(f)]);
disp(['TRIMGA violates ' num2str(sum(limits)) ' gene range limits']);
disp(' ');


function fitness = dummy(n,fixedIndex,variableIndex, ...
                         fixedValue,variableValue,P,GAP,D)
% DUMMY      an anonymous function that is used to call the given fitness 
%            function and execute the simplex method. This function is 
%            necessary in the case of two or more input arguments of the 
%            fitness function are required. In this case, it is required to 
%            avoid perturbation of the integer type parameters. 
%
% fitness  = dummy(n,fixedIndex,variableIndex,fixedValue, 
%                  variableValue,P,GAP,D)
%
% Inputs:
% n             = total number of design parameters
% fixedIndex    = the index of integer type parameters
% variableIndex = the index of the non-integer type parameters
% fixedValue    = the value of the integer parameters
% variableValue = the initial value for the non-integer parameters
% P             = Population
% GAP           = Genetic Algorithm Parameters
%                 (see gap_default for definitions)
% D             = data needed by fitness function (this could be an 
%                 empty matrix)
%
% Output:
% fitness       = fitness value
%
% Written by:
% Ricky Chan and Aaron Cramer for S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% get best individual
[~,bestFitIndex] = max(P.mfit(GAP.fp_obj,:));

% new vector containing both fixed and non-fixed parameters
value = nan(n,1);
value(fixedIndex) = fixedValue;
value(variableIndex) = variableValue;

% computing the fitness value
if nargin > 7 % if extra data is present
    % if extra information is passed to the fitness function evaluation
    if GAP.ev_ssd
        fitness = feval(P.fithandle,value,P.age(bestFitIndex),...
                        P.mfit(:,bestFitIndex),...
                        P.region(:,bestFitIndex),D);
    else 
        % if no additional information is 
        % required for the fitness function evaluation
        fitness = feval(P.fithandle,value,D);  
    end
else
    % if additional information is required for the 
    % fitness function evaluation
    if GAP.ev_ssd
        fitness = feval(P.fithandle,value,...
                        P.age(bestFitIndex),P.mfit(:,bestFitIndex),...
                        P.region(:,bestFitIndex));
    else
        fitness = feval(P.fithandle,value);
    end
end

fitness = fitness(GAP.fp_obj);

index=find((value'>GAP.gd_max)|(value'<GAP.gd_min));
if ~isempty(index)
   fitness=fitness-Inf;
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

