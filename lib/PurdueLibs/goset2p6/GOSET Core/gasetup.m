function [P,GAS] = gasetup(popsize,GAP,fithandle,D)
% GASETUP Sets up a population of chromosomes by defining initial values, 
%         and seting up and initializing the P data structure.  
%         This includes the initial evaluation of the fitness function
%
% [P,GAS] = gasetup(popsize,GAP,@fitness_function,[D])
%
% Inputs:
% popsize            = size (number of members of) the population
% GAP                = genetic algorithm parameters (see gap_default)
% fitness_function   = function to be used to evaluate the population
% D                  = an optional parameter containing information needed 
%                      by the fitness function
%
% Outputs:
% P   = Sturcture with initialize population of chromosomes
% GAS = Structure with genetic algorithm statistics
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

% set up fitness function
P.fithandle=fithandle;

% set the population size
if (round(popsize)==popsize)&&(popsize>0)
   P.size=popsize;
else
   error('Error in POPSETUP: Population size not proper');
end

% initialize fitness
P.mfit(1:GAP.fp_nobj,1:popsize)=NaN;
P.fit(1:popsize)=NaN;

% initialize evaluation status
P.eval(1:popsize)=0;

% initialize age of population
P.age(1:popsize)=0;

% check and initialize the fields
if min(size(GAP.gd_min)==size(GAP.gd_max))==0
   error(['Error in GASETUP:' ...
          ' minvals and maxvals not dimensionally compatible']);
end
if min(size(GAP.gd_min)==size(GAP.gd_type))==0
   error(['Error in GASETUP:' ...
          ' minvals and types not dimensionally compatible']);
end
if min(size(GAP.gd_min)==size(GAP.gd_cid))==0
   error(['Error in GASETUP:' ...
          ' minvals and chrom_id not dimensionally compatible']);
end
if size(GAP.gd_min,1)==1
   P.ngenes=size(GAP.gd_min,2);
else
   if size(GAP.gd_min,2)==1
      P.ngenes=size(GAP.gd_min,1);
   else
      error(['Error in GASETUP:' ...
             ' minvals and maxvals should be vectors; not arrays']);
   end
end
if min(GAP.gd_max>=GAP.gd_min)==0
   error(['Error in GASETUP:' ...
          ' element of minvals greater than' ...
          ' corresponding element of maxvals']);
end
P.min(1:P.ngenes,1)=GAP.gd_min(:);
P.max(1:P.ngenes,1)=GAP.gd_max(:);
P.type(1:P.ngenes,1)=GAP.gd_type(:);
P.chrom_id(1:P.ngenes,1)=GAP.gd_cid(:);
P = unrndinit(P,GAP);

% initialize the statistical information 
GAS.cg=1;
GAS.medianfit(1:GAP.fp_nobj,1:GAP.fp_ngen) = NaN;
GAS.meanfit(1:GAP.fp_nobj,1:GAP.fp_ngen)   = NaN;
GAS.bestfit(1:GAP.fp_nobj,1:GAP.fp_ngen)   = NaN;
GAS.bestgenes(1:P.ngenes,1:GAP.fp_ngen,1:GAP.fp_nobj) = NaN;
GAS.ne(1:GAP.fp_ngen) = NaN;

% initialize the fitness
P.mfit=NaN(GAP.fp_nobj,P.size);
if (nargin==3)
   [P.mfit,P.eval,GAS.ne(1)]=evaluate(P,GAP,0);
else
   if (nargin==4)
      [P.mfit,P.eval,GAS.ne(1)]=evaluate(P,GAP,0,D);
   else
      error('GASETUP: Incorrect number of input arguments');
   end
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
