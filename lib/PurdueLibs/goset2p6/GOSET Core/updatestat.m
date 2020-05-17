function [GAS] = updatestat(GAP,GAS,Pk)
% UPDATESTAT Update the statistic information of GAS
%
% [GAS] = updatestat(GAP,GAS,Pin)
%
% Inputs:
% GAP      = genetic algorithm parameter structure
% GAS      = a structure of genetic algorithm statistics
% Pin      = current population of chromosomes (a structure)
%
% Outputs:
% GAS      = a structure of genetic algorithm statistics
%
% Written by:
% Yonggon Lee
% United States Naval Academy
% E-mail: ylee@usna.edu

j = GAS.cg;
GAS.meanfit(:,j)  = mean(Pk.mfit,2);
GAS.medianfit(:,j) = median(Pk.mfit,2);
[GAS.bestfit(:,j), bestindex] = max(Pk.mfit,[],2);
GAS.bestgenes(1:Pk.ngenes,j,1:GAP.fp_nobj) = Pk.gene(:,bestindex');

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
