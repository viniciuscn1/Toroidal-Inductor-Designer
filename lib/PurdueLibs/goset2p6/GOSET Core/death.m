function [D_list] = death(P,P_list,GAP)
% DEATH   Creates a list of individuals to die (i.e. to be replaced by 
%         their children)
%
% [D_list] = death(P,P_list,GAP)
%
% Inputs:
% P        = current population of chromosomes (a structure)
% P_list   = parent list (see select)
% GAP      = genetic algorithm parameters
%
% Outputs:
% D_list  =  Death list.   This is an array wherein the rows correspond 
%            to different regions and whose columns are the indcices of the 
%            individuals in that retion to die. Note, values with zero 
%            indicate that the corresponding region has fewer individuals
%            than the region with the max (the zeros are used for packing 
%            purposes)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize death list
D_list=zeros(size(P_list));

% determine the number of regions
N_regions=size(D_list,1);

% determine number of parents in each region
for region=1:N_regions,
   index=find(P_list(region,:)==0,1,'first');
   if ~isempty(index)
      D_size(region)=index-1;
   else
      D_size(region)=size(D_list,2);
   end
end

% now determine the death list
if GAP.dt_alg==6
   temp=randperm(4);
   death_algorithm=temp(1);
else
   death_algorithm=GAP.dt_alg;
end

switch death_algorithm
       
   case 1,             % replace parents with children
       D_list=P_list;
           
   case 2,             % randomly select death list
      for region=1:N_regions,
         region_list=find(P.region==region);
         region_size=length(region_list);
         randomized_region_list=region_list(randperm(region_size));
         D_list(region,1:D_size(region))= ...
               randomized_region_list(1:D_size(region));
      end
      
    case 3,            % death by tournament selection on aggregate fitness
       for region=1:N_regions,
          region_list=find(P.region==region);
          region_size=length(region_list);
          region_fit=P.fit(region_list);
          for i=1:D_size(region),
             mworst=ceil(max([eps rand])*region_size);
             mworstfit=region_fit(mworst);
             for n=2:GAP.dt_nts,
                m=ceil(max([eps rand])*region_size);
                if region_fit(m)<mworstfit
                   mworst=m;
                   mworstfit=region_fit(m);
                end
             end
             D_list(region,i)=region_list(mworst);
          end
       end
       
    case 4,           % death by tournament selection on age (oldest dies)
       for region=1:N_regions,
          region_list=find(P.region==region);
          region_size=length(region_list);
          region_age=P.age(region_list);
          for i=1:D_size(region),
             moldest=ceil(max([eps rand])*region_size);
             moldestage=region_age(moldest);
             for n=2:GAP.dt_nts,
                m=ceil(max([eps rand])*region_size);
                if region_age(m)>moldestage
                   moldest=m;
                   moldestage=region_age(m);
                end
             end
             D_list(region,i)=region_list(moldest);
          end
       end
      
    case 5,            % custom algorithm
       for region=1:N_regions,
          region_list=find(P.region==region);
          region_age=P.age(region_list);
          region_mfit=P.mfit(region_list);
          region_fit=P.fit(region_list);
          d_list=feval(GAP.dt_cah,region,D_size(region), ...
                       region_age,region_mfit,region_fit);
          D_list(region,:)=region_list(d_list);
       end 
      
      
   otherwise
      
      error('GAP.dt_alg is not valid');
          
end % switch
           
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