function [O] = mutate(N,GAP)
% MUTATE  Perform genetic mutation on a population of chromosomes.
%         Four mutation algorithms are sequentially applied
%         1.) Total mutationa -genes that are mutated have no
%             relationship to their previous value
%         2.) Partial mutation - all non-integer genes that are mutated 
%             are perturbed from their original value (both relative and 
%             absolute mutations are introduced)
%         3.) Vector mutation - all non-interger genes on a chromosome are
%             perturbed in a random direction (both relative and absolute 
%             mutations are introduced)
%         4.) Integer mutation - an extra total mutation operator for 
%             integer genes.
%
% [O] = mutate(N,GAP)
%
% Inputs:
% N        = current population of chromosomes (a structure)
% GAP      = structure of genetic algorithm parameters
%
% Outputs:
% O        = population of chromosomes after mutation (a structure)
%
% Written by:
% S.D. Sudhoff
% Scott Sudhoff
% Purdue University
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize mutated set
O=N;

% Total Mutation ---------------------------------------------------------%
if GAP.mt_ptgm>0

   % total number of genes to be mutated
   Nmutations=round(GAP.mt_ptgm*O.ngenes*O.size);
   if (Nmutations < 1)
      Nmutations=1;
   end

   for i=1:Nmutations,
    
      % identify member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % identify gene to be mutated
      gene=round(rand*(O.ngenes-1)+1);
      
      % mutate the gene
      if O.type(gene)==1
         levels=N.max(gene)-N.min(gene)+1;
         if levels>1
            O.normgene(gene,member)=(fix(rand*levels))/ ...
                                    (levels-1);
         else
            O.normgene(gene,member)=0;
         end
      else                       
         O.normgene(gene,member)=rand;
      end % switch
         
      % mark that it must be reevaluated
      O.eval(member)=0;
          
   end % mutation loop
   
end % if


% Partial Mutation -------------------------------------------------------%

% Relative Partial Mutation
if (GAP.mt_prgm>0)&&(GAP.mt_srgm>0)

   % total number of genes to be mutated
   Nmutations=round(GAP.mt_prgm*O.ngenes*O.size);
   if (Nmutations < 1)
      Nmutations=1;
   end

   for i=1:Nmutations,
    
      % identify population member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % identify gene to be mutated
      gene=floor(rand*O.ngenes+1.0-1.0e-12);
      
      % mutate the gene
      if O.type(gene)>1
         O.normgene(gene,member)=O.normgene(gene,member)* ...
                                 (1+randn*GAP.mt_srgm);
      end
       
      % limit the range
      O.normgene(gene,member)=generepair(O.normgene(gene,member),GAP);
            
      % mark that it must be reevaluated
      O.eval(member)=0;
      
   end % for loop on mutations
    
end

% Absolute Partial Mutation
if (GAP.mt_pagm>0)&&(GAP.mt_sagm>0)

   % total number of genes to be mutated
   Nmutations=round(GAP.mt_pagm*O.ngenes*O.size);
   if (Nmutations < 1)
      Nmutations=1;
   end

   for i=1:Nmutations,
    
      % identify population member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % identify gene to be mutated
      gene=floor(rand*O.ngenes+1.0-1.0e-12);
      
      % mutate the gene
      if O.type(gene)>1
         O.normgene(gene,member)=O.normgene(gene,member)+randn*GAP.mt_sagm;
      end
       
      % limit the range
      O.normgene(gene,member)=generepair(O.normgene(gene,member),GAP);
           
      % mark that it must be reevaluated
      O.eval(member)=0;
       
   end % for loop on mutations
   
end

% Vector Mutation --------------------------------------------------------%

% Relative Vector Mutation
if (GAP.mt_prvm>0)&&(GAP.mt_srvm>0)&&(min(N.type)>1)

   % total number of genes to be mutated
   Nmutations=round(GAP.mt_prvm*O.size);
   if (Nmutations < 1)
      Nmutations=1;
   end

   for i=1:Nmutations,
    
      % identify population member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % mutate the gene
      direction=randn([O.ngenes 1]).*(O.type>1);
      ndirection=direction/norm(direction);
      O.normgene(:,member)=O.normgene(:,member).* ...
                           (1+ndirection*GAP.mt_srvm*randn);
       
      % limit the range
      O.normgene(:,member)=generepair(O.normgene(:,member),GAP);
     
      % mark that it must be reevaluated
      O.eval(member)=0;
      
   end % for loop on mutations
   
end

% Absolute Vector Mutation
if (GAP.mt_pavm>0)&&(GAP.mt_savm>0)&&(min(N.type)>1)

   % total number of genes to be mutated
   Nmutations=round(GAP.mt_pavm*O.size);
   if (Nmutations < 1)
      Nmutations=1;
   end

   for i=1:Nmutations,
    
      % identify population member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % mutate the gene
      direction=randn([O.ngenes 1]).*(O.type>1);
      ndirection=direction/norm(direction);
      O.normgene(:,member)=O.normgene(:,member)+ ...
                           ndirection*GAP.mt_savm*randn;
       
      % limit the range
      O.normgene(:,member)=generepair(O.normgene(:,member),GAP);
      
      % mark that it must be reevaluated
      O.eval(member)=0;
      
   end % for loop on mutations
    
end

% Integer Mutation -------------------------------------------------------%

% find list of integer genes
iglist=find(O.type==1);

% find number of integer genes
niglist=length(iglist);

% find number of mutations
Nmutations=round(GAP.mt_pigm*niglist*O.size);

% decide if mutation needed
if Nmutations>0

   for i=1:Nmutations,
    
      % identify member to be mutated
      member=ceil(max([eps rand])*O.size);
   
      % identify gene to be mutated
      gene=iglist(round(rand*(niglist-1)+1));
           
      % mutate the gene
      levels=N.max(gene)-N.min(gene)+1;
      if levels>1
         O.normgene(gene,member)=(fix(rand*levels))/(levels-1);
      else
         O.normgene(gene,member)=0;
      end
      
      % mark that it must be reevaluated
      O.eval(member)=0;
          
   end % mutation loop
    
end % if

% Now that we are done mutating the heck out of everything, 
% lets compute the raw gene values
O=rawgene(O);

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
