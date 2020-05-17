function [P_list,PL_size] = select(P,GAP)
% SELECT  Selects a mating pool from a population of chromosomes.
%         Note that the most fit chromosome is always selected if elitism
%         is used.
%
% [P_list,PL_size] = select(P,GAP)
%
% Inputs:
% P        = current population of chromosomes (a structure)
% GAP      = genetic algorithm parameters
%
% Outputs:
% P_list  =  Parent list. This is an array wherein the rows correspond to
%            different regions and whose columns are the indcices of the 
%            individuals in the mating pool. Note, values with zero 
%            indicate that the corresponding region has fewer individuals
%            than the region with the max (the zeros are used for packing
%            purposes)
% PL_size =  Parent list size.  Elements describe the number of parents in 
%            each region.
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

% check to make sure all chromosomes have been evaluated
if min(P.eval==0)
   error('Error in Select: Not all chromosomes have been evaluated');
end

% find the number of regions
N_regions=max(P.region);

% initialize parent list size
PL_size=zeros(N_regions,1);

% determine number of parents in each region
for region=1:N_regions,
   region_list=find(P.region==region);
   region_size=length(region_list);
   PL_size(region)=round(region_size*GAP.mc_pp/2)*2;
   if PL_size(region)<2
      PL_size(region)=2;
   end
end

% initialize parent list
P_list(N_regions,max(PL_size))=0;

% perform selection separately for each region
for region=1:N_regions,
       
    % find the regional population and its fitness after niching
    region_list=find(P.region==region);
    region_fit=P.fit(region_list);
   
    switch GAP.sl_alg
        
       case 1 % roulette wheel selection
          % determine the mating probability
          matprob=region_fit/sum(region_fit);
          % create a mapping function for selection
          map=cumsum(matprob);
          % now do the selection
          for i=1:PL_size(region),
             choice=rand;
             j=1;
             while (choice > map(j))
                j=j+1;
             end
             P_list(region,i)=region_list(j);
          end
            
       case 2 % tournament selection
          for i=1:PL_size(region),
             mbest=ceil(max([eps rand])*PL_size(region));
             mbestfit=region_fit(mbest);
             for n=2:GAP.sl_nts,
                m=ceil(max([eps rand])*PL_size(region));
                if region_fit(m)>mbestfit
                   mbest=m;
                   mbestfit=region_fit(m);
                end
             end
             P_list(region,i)=region_list(mbest);
          end 
          
       case 3 % custom selection
          region_mfit=P.mfit(region_list);
          region_age=P.age(region_list);
          p_list=feval(GAP.sl_cah,region,PL_size(region), ...
                       region_age,region_mfit,region_fit);
          P_list(region,:)=region_list(p_list);
          
       otherwise
          error('Invalid selection option');
                 
    end % switch
       
end  % over each region 

% References
% [1] Genetic Algorithms in Search, Optimization, and Machine Learning
% by D.E. Goldberg, Adison Wesley Longman, 1989

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
