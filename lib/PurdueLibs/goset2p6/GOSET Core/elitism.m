function [O] = elitism(N1,N0,GAP,GAS)
% ELITISM performs an elitism operation by ensuring that a population
%         that has been processed by several genetic operators will
%         have an individual in each region that is at least fit as the 
%         fittest individual in each region prior to the application of
%         the genetic operators.  In multi-objective optimizations this
%         routine will ensure that no non-dominated members of the 
%         population will be lost, if possible.
%
% [O] = elitism(N1,N0,GAP,GAS)
%
% Inputs:
% N1       = processed population of chromosomes (a structure)
% N0       = original population of chromosomes (a structure)
% GAP      = structure of genetic algorithm parameters
% GAS      = genetic algorithm statistics
%
% Outputs:
% O        = population of chromosomes after elitism (a structure)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize the output to the new population
O=N1;


% multi-objective elitism
if (GAP.el_act==1)&&(GAP.fp_nobj>1)&& ...
   (GAP.fp_obj==0)&&(GAS.cg>GAP.el_fgs*GAP.fp_ngen)

   % perform elitism on each region
   for i=1:max(N1.region)
       
      % or_list is the list of individuals in 
      % the old population in the region of interest
      or_list=find(N0.region==i);
      
      % remove redundant individuals from old region list
      [~,imap,~]=unique(N0.mfit(:,or_list)','rows');
      or_list=or_list(imap);
      
      % remove dominated individuals from old region list
      nd=nondom(N0.mfit(:,or_list),1);
      or_list=or_list(nd==1);
      
      % nr_list is the list of individuals in the 
      % new population in the region of interest
      nr_list=find(N1.region==i);
      nr_list2=nr_list;
      
      % compute number of individuals in region
      nr=length(nr_list);
        
      % remove redundant individuals from new region list
      [~,imap,~]=unique(N1.mfit(:,nr_list)','rows');
      nr_list=nr_list(imap);
            
      % undnr_list is the the list of unique, non-dominated 
      % individuals in the new population in this region
      nd=nondom(N1.mfit(:,nr_list),1);
      nr_list=nr_list(nd==1);
      
      % new remove members from old regional 
      % list who are in the new region list
      [~,map]=setdiff(N0.mfit(:,or_list)', N1.mfit(:,nr_list)','rows');
      or_list=or_list(map);
            
      % hp is a hybrid population of unique, non-dominated individuals 
      % of the old and new populations
      hp_mfit =[N0.mfit(:,or_list) N1.mfit(:,nr_list)];

      % see which elements of hp are non-dominated
      hp_nd=nondom(hp_mfit,1);
      
      % modify the or_list and nr_list to include only individuals who are 
      % non-dominated in hybrid population
      temp1=length(or_list);
      temp2=length(hp_nd);
      or_list=or_list( hp_nd(1:temp1)==1          );
      nr_list=nr_list( hp_nd((temp1+1):temp2)==1  );
      
      % compute final number of individuals in old population 
      % who will be non-dominated in hybrid population
      no=length(or_list);
      
      % compute final number of individuals in new population who 
      % will be non-dominated in hybrid population
      nn=length(nr_list);
     
      % compute number of non-dominated individuals in hybrid population
      nh=no+nn;
      
      % compute number of non-dominated solutions to be kept
      nk=round(nr*GAP.el_fpe);
      
      if (nh<nk)
         
         % put all members of or_list into new population
         % delete dominated members of new population at random
         % to make room
         save_list=or_list;
         delete_list=setdiff(nr_list2,nr_list);
         delete_list=delete_list(randperm(length(delete_list)));
         delete_list=delete_list(1:length(save_list));
                   
      else
         
         % keep those members of hybrid population with best penalty
         % compute penalty associated with hybrid population
         hp_mfit=hp_mfit(:,hp_nd==1);
         hp_penalty=crowding(hp_mfit);
         [~,index]=sort(-hp_penalty);
         
         % list of the solutions we will presearve as indexed in hp_sfit
         preserve_list=index(1:nk);
         
         % find list of solutions to save from the old population
         save_list=or_list(index(preserve_list<=no));
         
         % find list of solutions to save from new population
         % these are saved by virtue of being exempt from the delete list
         save_list2=nr_list(index(preserve_list>no)-no);
         
         % now form the delete list
         delete_list=setdiff(nr_list2,save_list2);
         delete_list=delete_list(randperm(length(delete_list)));
         delete_list=delete_list(1:length(save_list));
       
      end
 
      % now substitute until we run out of old non-dominated solutions
      % or places to put them within the region.  
      O.gene(:,delete_list)     = N0.gene(:,save_list);
      O.normgene(:,delete_list) = N0.normgene(:,save_list);
      O.mfit(:,delete_list)     = N0.mfit(:,save_list);
      O.age(delete_list)        = N0.age(save_list)+1;
      
   end  % for loop over regions
   
end % multi-objective elitism

% single objective elitism
% if the older population has a individual better than the new population,
% replace the worst individual in the new population with the best 
% iindividual in the old population
if (GAP.el_act==1)&&(GAS.cg>GAP.el_fgs*GAP.fp_ngen)
   
   for i=1:max(N1.region)
       
      % or_list is the list of individual in a 
      % given region in the old population
      or_list=find(N0.region==i);
      
      % nr_list is the list of individuals in the 
      % given region in the new population
      nr_list=find(O.region==i);
      
      for o=1:GAP.fp_nobj
      
         % or_best,or_best_index are the best fitness and index of 
         % individual with best fitness in or_list
         [or_best,or_best_index]=max(N0.mfit(o,or_list));
            
         % nr_best and nr_best_index are the best fitness and index of 
         % individual with the best index in nr_list
         [nr_best,~]=max(O.mfit(o,nr_list));
      
         % nr_worst and nr_worst_index are the worst fitness and index of 
         % individual with the worst index in nr_list
         [~,nr_worst_index]=min(O.mfit(o,nr_list));
      
         % replace as necessary
         if (or_best>nr_best)
            oi=or_list(or_best_index);
            ni=nr_list(nr_worst_index);
            O.gene(:,ni)=N0.gene(:,oi);
            O.normgene(:,ni)=N0.normgene(:,oi);
            O.mfit(:,ni)=N0.mfit(:,oi);
            O.age(ni)=N0.age(oi)+1;
         end
        
         
      end % for over number of objectives  
      
   end % for over number of regions
   
end % single objective elitism


function p=crowding(x)

   % determine the number of vectors n, 
   % and number of elements in each vector, e
   [e,n]=size(x);
 
   for i=1:e,
      emax=max(x(i,:));
      emin=min(x(i,:));
      if (emax~=emin)
         x(i,:)=(x(i,:)-emin)/(emax-emin);
      else
         x(i,:)=0.5;
      end
   end
   
   % initialize the penalty
   p=ones(1,n);
      
   % compute the distance between every element and every other element
   d(n,n)=0.0;
   for i=1:n,
      d(i,i)=sqrt(e);
      xi=x(:,i);
      for j=i+1:n,
         xj=x(:,j);
         d(i,j)=norm(xj-xi);
         d(j,i)=d(i,j);
      end
   end
      
   % determine the penalty
   for i=1:n,
       p(i)=min(d(i,:));
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