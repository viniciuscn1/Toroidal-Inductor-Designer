function [fit] = scale(P,GAP)
% SCALE      Given the current genetic algorithm parameter set (GAP)
%            and a current population (P) this routine computes the
%            scaled and aggregated fitness.
%
% [fit] = scale(P,GAP)
%
% Inputs:
% P        = current population of chromosomes (a structure)
% GAP      = current value of genetic algorithm parameters
%
% Outputs:
% fit      = scaled and aggregated fitness
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

% initialize
fit(1:P.size)=0;

% check to make sure all chromosomes have been evaluated
if min(P.eval==0)
   error('Error in Select: Not all chromosomes have been evaluated');
end

% find the number of objectives
nobj=size(P.mfit,1);

% find the number of regions
nreg=max(P.region);

% perform scaling separately for each region
for region=1:nreg,
 
    % find the regional population
    region_list=find(P.region==region);
    
    % declare scaled fitness
    scaledfitness(nobj,length(region_list))=0;
    
    % now scale each objective
    for obj=1:nobj
    
       % compute the fitness of the individuals in the region
       region_fit=P.mfit(obj,region_list);
              
       % look at fitness within region
       fmin=min(region_fit);
       fmax=max(region_fit);
       
       if fmax>fmin       % case where not all fitness values are the same

          %adjust the fitness
          a=0.0;
          b=0.0;
          c=0.0;
          d=0.0; 
          switch(GAP.sc_alg)
             case{0} % no scaling
                a=1.0;
                b=0.0;
             case{1} % offset scaling
                a=1;
                b=-fmin;
             case{2} % standard linear scaling [1]
                % in this algorithm, the scaled fitness
                % has the same average as the original fitness
                % the maximum fitness is GAP.sc_kln more likely
                % to be selected than the average 
                % if all values are positive
                favg=mean(region_fit);
                a=(GAP.sc_kln-1.0)*favg/(fmax-favg);
                b=favg*(1-a);
              case{3} % modified linear scaling 
                % same as 2 above, but median is used rather than average
                % modified on 5/5/2011
                fmed=median(region_fit);
                a=(GAP.sc_kln-1.0)/(fmax-fmed);
                b=(1-a*fmed);
             case{4} % mapped linear scaling
                a=(GAP.sc_kln-1.0)/(fmax-fmin);
                b=-fmin*a+1.0;
             case{5} % sigma truncation of [1]
                favg=mean(region_fit);
                fstd=std(region_fit);
                a=1.0;
                b=-(favg-GAP.sc_cst*fstd);
             case{6} % quadratic scaling
                favg=mean(region_fit);
                A=[fmin*fmin fmin 1; favg*favg favg 1; fmax*fmax fmax 1];
                B=[GAP.sc_kmnq; 1.0; GAP.sc_kmxq];
                PV=A\B;
                a=PV(1);
                b=PV(2);
                c=PV(3);
          end % switch statement
       
          % first step in computing the scaled fitness 
          switch(GAP.sc_alg)
             case{0,1,2,3,4,5} % linear scaling
                newfitness=a*region_fit+b;
             case{6} % quadratic scaling
                newfitness=a*region_fit.*region_fit+b*region_fit+c;   
          end % switch statement
        
          % now clip negative values
          newfitness(newfitness<0)=0;
       
          % find the sum of new fitness
          sumnewfitness=sum(newfitness);
       
          % now assign the final scaled fitness, which would be the 
          % probability of an unpenalized (with respect to diversity 
          % control) single objective optimization
          scaledfitness(obj,region_list)=newfitness/sumnewfitness;
                            
       else

          % cannot really scale since all function values are the same
          scaledfitness(obj,region_list)=1.0/length(region_list);         
           
       end
       
   end % over each objective   
   
   % now also find the aggregate fitness
   fit(region_list)=(GAP.owv*scaledfitness(:,region_list)).* ...
                    (P.pen(region_list));  

end  % over each region 

% References
% [1] Genetic Algorithms in Search, Optimization, and Machine Learning
%     by D.E. Goldberg, Adison Wesley Longman, 1989

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