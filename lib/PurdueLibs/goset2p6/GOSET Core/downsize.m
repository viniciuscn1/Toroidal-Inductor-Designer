function [M] = downsize(P,newsize)
% DOWNSIZE  Creates a population of individuals based on the most fit
%           individuals from P.
%
% [M] = downsize(P,newsize)
%
% Inputs:
% P        = current population of chromosomes (a structure)
% newsize  = the size of the new population
%
% Outputs:
% M        = new populsation of chromosomes
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

% check on population size
if newsize>=P.size
   error('Error in GA Toolbox downsize: requested size >= current size');
end

% make sure the entire population has been evaluated
if min(P.eval==0)
   error('Error in GA Toolbox downsize: population not fully evaluated');
end

% copy the common fields
M.ngenes=P.ngenes;
M.size=newsize;
M.fithandle=P.fithandle;
M.eval(1:newsize)=1;
M.chrom_id=P.chrom_id;
M.min=P.min;
M.max=P.max;
M.type=P.type;

% number of population members left to assign and left to consider
nleft_new=newsize;
nleft_old=P.size;

% index of starting point to be assigned
j=1;

% assign region by region
for region=1:max(P.region),
         
   % compute a list of population members in the region 
   clear region_list;
   region_list=find([P.region==region]);
  

      % determine numer of populations member in the region
      % and number of population member in new population in this region
      n_old=length(region_list);
      n_new=round(nleft_new*(n_old/nleft_old));
      nleft_new=nleft_new-n_new;
      nleft_old=nleft_old-n_old;
   
      % sort the list of values by cumulative rank
      % the cumulative rank is the sum of the ranks over each objectives.
      % the best rank possible is obtained if the individual was the best
      % in very single objective
      clear rank;
      rank(1:n_old)=0;
      for obj=1:size(P.mfit,1),
        [~,obj_indices]=sort(P.mfit(obj,region_list));
        for k=1:n_old,
           rank(obj_indices(k))=rank(obj_indices(k))+k;
        end
      end
      [~,indices]=sort(rank);
    
      % compute locations
      old_elements=region_list(indices(n_old-n_new+1:n_old));
      new_elements=j:j+n_new-1;
      j=j+n_new;
   
      % assign elements as appropriate
      M.mfit(:,new_elements)     = P.mfit(:,old_elements);
      M.fit(:,new_elements)      = P.fit(:,old_elements);
      M.gene(:,new_elements)     = P.gene(:,old_elements);
      M.normgene(:,new_elements) = P.normgene(:,old_elements);
      M.region(new_elements)     = P.region(old_elements);
      M.age(new_elements)        = P.age(old_elements);         
      
end % region loop

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

