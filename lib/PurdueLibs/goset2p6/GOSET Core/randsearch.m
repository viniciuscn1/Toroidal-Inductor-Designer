function [O,GAS] = randsearch(N,GAP,GAS,D)
% Randsearch  This routine performs a random search in the vicinity of the
%             most fit individual in each region.
%
% [O,GAS] = randsearch(N,GAP,GAS,D)
%
% Inputs:
% N        = current population of chromosomes (a structure)
% GAP      = genetic algorithm parameters (see gap_default)
% GAS      = Genetic Algorithm Statistics
% D        = an optional data structure if needed for fitness evaluation
%
% Outputs:
% O        = new population of individuals with most fit individuals
%            potentially replaced by randsom search results
% GAS      = updated genetic algorithm statistics
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize the output
O=N;

% decide if we need to do a random search
if ((GAP.rs_srp>0)||(GAP.rs_sap>0))&&(GAP.rs_fps>0)&& ...
   (GAS.cg>GAP.rs_fgs*GAP.fp_ngen)&&(GAP.fp_obj>0)&&(GAP.rs_fea>rand)&& ...
   (max(N.type)>1)

   % determine the number of regions 
   nregions=max(N.region); 
   
   % determine the number of individuals to come out of each region
   % and the number of individuals used for the random search from that
   % region
   regional_population(1:nregions)=0;
   regional_nsearch(1:nregions)=0;
   regional_bestfit(1:nregions)=0;
   regional_bestfit_index(1:nregions)=0;
   for i=1:nregions
       region_list=find(N.region==i);
       regional_population(i)=max(size(region_list));
       [regional_bestfit(i),regional_bestfit_subindex]= ...
                                      max(N.mfit(GAP.fp_obj,region_list));
       regional_bestfit_index(i)=region_list(regional_bestfit_subindex);
       regional_nsearch(i)=ceil(GAP.rs_fps*regional_population(i));
   end
   
   % total number of search points
   total_nsearch=sum(regional_nsearch); 

   % Population to be used in search
   S.ngenes=N.ngenes;
   S.gene(1:N.ngenes,1:total_nsearch)=0;
   S.normgene(1:N.ngenes,1:total_nsearch)=0;
   S.min=N.min;
   S.max=N.max;
   S.eval(1:total_nsearch)=0;
   S.fit(1:total_nsearch)=NaN;
   S.mfit(1:GAP.fp_nobj,1:total_nsearch)=NaN;
   S.region(1:total_nsearch)=1;
   S.fithandle=N.fithandle;
   S.size=total_nsearch;
   S.type=N.type;
     
   if GAP.rs_frp>rand
  
      % conduct randsom search with relative perturbations
      l=1;
      for i=1:nregions
         % identify the best individual
         bestgenes=N.gene(:,regional_bestfit_index(i));
         % create a list of genes
         for j=1:regional_nsearch(i)
            radious=randn*GAP.rs_srp;
            direction=randn(size(bestgenes)).*(N.type>1);
            ndirection=direction/norm(direction);
            S.gene(:,l)=bestgenes.*(1+ndirection*radious);
            S.age(l)=N.age(regional_bestfit_index(i));
            for k=1:max(size(bestgenes))
               if S.gene(k,l) > N.max(k)
                  S.gene(k,l) = N.max(k);
               end
               if S.gene(k,l) < N.min(k)
                 S.gene(k,l) = N.min(k);
               end
            end
            S.region(l)=i;
            l=l+1;
         end
      end % loop over regions
      % update the actual value
      S=normgene(S);
      % have completed individuals with relative perturbations
      
   else
 
      % conduct randsom search with absolute perturbations
      l=1;
      for i=1:nregions
         % identify the best individual
         bestgenes=N.normgene(:,regional_bestfit_index(i));
         % create a list of genes
         for j=1:regional_nsearch(i)
            radious=randn*GAP.rs_sap;
            direction=randn(size(bestgenes)).*(N.type>1);
            ndirection=direction/norm(direction);
            S.normgene(:,l)=bestgenes+ndirection*radious;
            S.age(l)=N.age(regional_bestfit_index(i));
            for k=1:max(size(bestgenes))
               if S.normgene(k,l) > 1
                  S.normgene(k,l) = 1;
               end
               if S.normgene(k,l) < 0
                  S.normgene(k,l) = 0;
               end
            end
            S.region(l)=i;
            l=l+1;
         end
      end % loop over regions
      % update the actual value
      S=rawgene(S);
      % have completed individuals with absolute perturbations       
           
   end
   
   % do an evaluation (and scaling and aggregation) 
   % based on the random population
   [S.mfit,S.eval,GAS.ne(GAS.cg)]=evaluate(S,GAP,GAS.ne(GAS.cg),D);
   
   % if we found any better solutions, put them into the population
   for i=1:nregions
       regional_search_list=find(S.region==i);
       [regional_search_max,regional_search_sub_index]= ...
                             max(S.mfit(GAP.fp_obj,regional_search_list));
       regional_search_index= ...
                           regional_search_list(regional_search_sub_index);
       if regional_search_max>regional_bestfit(i)
          O.gene(:,regional_bestfit_index(i))= ...
                                           S.gene(:,regional_search_index);
          O.normgene(:,regional_bestfit_index(i))= ...
                                       S.normgene(:,regional_search_index);
          O.mfit(:,regional_bestfit_index(i))= ...
                                           S.mfit(:,regional_search_index);
          O.age(regional_bestfit_index(i))=0;
      end
   end
   
end % conduct the random search

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