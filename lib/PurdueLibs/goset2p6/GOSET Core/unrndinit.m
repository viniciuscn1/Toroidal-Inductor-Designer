function [P] = unrndinit(P,GAP)
% UNRNDINIT Uniform random initialization of a gene pool. Initializes
%           the value and region of the genes
%
% [P] = unrndinit(P,GAP)
%
% Inputs:
% P        = current population of chromosomes (a structure)
% GAP      = a structure with the genetic algorithm parameters (see
%            gap_default)
%
% Outputs:
% P        = population of chromosomes after genes randomly 
%            initialized (a structure)
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

% initialize the values of the genes
P.normgene(P.ngenes,P.size)=0;
for j=1:P.ngenes,
   if P.type(j)==1
      levels=P.max(j)-P.min(j)+1;
      if levels>1
         P.normgene(j,:)=(fix(rand(1,P.size)*levels))/ ...
                         (levels-1);
      else
         P.normgene(j,:)=0;
      end
   else
      P.normgene(j,:)=rand(1,P.size); 
   end
end

% compute the raw values of the gene
P=rawgene(P);

% distribute into regions
P.region=ceil(max(eps,rand(1,P.size))*GAP.mg_nreg);

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
