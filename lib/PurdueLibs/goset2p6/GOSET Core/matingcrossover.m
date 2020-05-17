function [N] = matingcrossover(P,P_list,PL_size,D_list,GAP,GAS)
% MATINGCROSSOVER Perform mating and genetic crossover on a population.
%                 The most fit chromosome is not allowed to be replaced
%
% [N] = matingcrossover(P,P_list,PL_size,D_list,GAP,GAS)
%
% Inputs:
% P        = current population (a structure)
% P_list   = parent list (see select)
% PL_size  = parent list size
% D_list   = death list (see death)
% GAP      = a structure containing the parameters of the genetic algorithm
%            (see gap_default)
% GAS      = Genetic Algorithm Statistics
%
% Outputs:
% N        = population of chromosomes after chrossover (a structure)
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

% initialize output population
N=P;

% compute number of regions
N_regions=max(P.region);

% determine the crossover alogrithm 
persistent ca;
if GAP.mc_alg~=6 
   ca=GAP.mc_alg;
else   
   if (GAS.cg==2)
      ca=-1; 
   else
      if (GAS.cg>GAP.mc_gac+1)
         if GAS.bestfit(GAS.cg-1)==GAS.bestfit(GAS.cg-1-GAP.mc_gac)
            ca=-1;
         end
      end    
   end
   if isempty(ca)
      ca=-1;
   end
   if (ca==-1)
      ca=round(0.5+rand*5);
      if ca<1
         ca=1;
      end
      if ca>5
         ca=5;
      end
   end
end

% perform crossovers for each population subset
for region=1:N_regions,
 
   % find the number children needed in the region
   N_children=PL_size(region);
   

   for j=1:ceil(N_children/2),
       
      % determine the parent and children indices
      j1=2*j-1;
      j2=j1+1;
      child1=D_list(region,j1);
      child2=D_list(region,j2);
      parent1=P_list(region,j1);
      parent2=P_list(region,j2);

      % null out the members that will be crossed over
      N.eval(child1)=0;
      N.eval(child2)=0;
     
      % exchange genetic material
      switch ca 
         case {1}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            single_point_crossover(P.normgene(:,parent1), ...
                                   P.normgene(:,parent2), ...
                                   P.chrom_id,GAP.mc_fc);
         case {2}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            scalar_simple_blend_crossover(P.normgene(:,parent1), ...
                                          P.normgene(:,parent2), ...
                                          P.chrom_id,P.type,GAP);
         case {3}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
             vector_simple_blend_crossover(P.normgene(:,parent1), ...
                                           P.normgene(:,parent2), ...
                                           P.chrom_id,P.type,GAP);
         case {4}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            scalar_simulated_binary_crossover(P.normgene(:,parent1), ...
                                              P.normgene(:,parent2), ...
                                              P.chrom_id,P.type,GAP);
         case {5}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            vector_simulated_binary_crossover(P.normgene(:,parent1), ...
                                              P.normgene(:,parent2), ...
                                              P.chrom_id,P.type,GAP);
         case {7}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            multi_point_crossover(P.normgene(:,parent1), ...
                                  P.normgene(:,parent2), ...
                                  P.chrom_id,0.33);
         case {8}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            mpsbc(P.normgene(:,parent1),P.normgene(:,parent2), ...
                  P.chrom_id,GAP,0.33);       
        
         case {9}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            scalar_quasi_biological_crossover(P.normgene(:,parent1), ...
                                              P.normgene(:,parent2), ...
                                              P.chrom_id,P.type,GAP);
         case {10}
            [N.normgene(:,child1),N.normgene(:,child2)]= ...
            vector_quasi_biological_crossover(P.normgene(:,parent1), ...
                                              P.normgene(:,parent2), ...
                                              P.chrom_id,P.type,GAP);
         otherwise
            error('invalid crossover routine'); 
      end
         
   end % loop over number of individuals
   
end % loop over number of regions

% update raw gene values
N=rawgene(N);

% single point crossover based mating
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% fc       =  fraction of chromosomes to be crossed over
% c1       =  genes of child 1
% c2       =  genes of child 2
function [c1,c2]=single_point_crossover(p1,p2,chrom_id,fc)

   % initialize the children
   c1=p1;
   c2=p2;
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (fc>rand)&&(N_genes>1)
            
         % determine the location of the crossover
         location=ceil(max([eps rand])*(N_genes-1));
   
         % get the genetic material for these chromosomes
         range=1:location; 
         c1(gene_list(range))=p2(gene_list(range));
         c2(gene_list(range))=p1(gene_list(range));
               
      else % randomly swap the entire chromosome
  
         % assign gene values (different cromosomes are segregated 
         % but not crossed over) 
         if rand>0.5
            c1(k)=p2(k);
            c2(k)=p1(k);
         end   
          
      end % crossover versus swap
        
  end % loop over the number of chromosomes
  
% multi point crossover based mating
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% fc       =  probability of crossover at gene boundary
% c1       =  genes of child 1
% c2       =  genes of child 2
function [c1,c2]=multi_point_crossover(p1,p2,chrom_id,pc)

   % initialize the children
   c1=p1;
   c2=p2;   

   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));

      % initialize children
      swap=2.0*(rand>0.5)-1;
      gene=gene_list(1);
      if swap>0
         c1(gene)=p1(gene);
         c2(gene)=p2(gene);
      else
         c1(gene)=p2(gene);
         c2(gene)=p1(gene);
      end
      
      % due the crossover
      for l=1:N_genes-1,
         gene=gene_list(l);
         a=0.5*(1.0+swap);
         b=0.5*(1.0-swap);
         c1(gene)= p1(gene)*a+p2(gene)*b;
         c2(gene)= p2(gene)*a+p1(gene)*b;
         if pc>rand
            swap=-swap;
         end
      end
      
end % loop over the number of chromosomes

% multi point simulated binary crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% fc       =  probability of crossover at gene boundary
% c1       =  genes of child 1
% c2       =  genes of child 2
function [c1,c2]=mpsbc(p1,p2,chrom_id,GAP,pc)

   % initialize the children
   c1=p1;
   c2=p2;   

   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));

      % initialize the swap
      swap=2.0*(rand>0.5)-1;
      a=0.5*(1+swap);
      b=0.5*(1-swap);
      
      % due the crossover
      for l=1:N_genes,
         gene=gene_list(l);
         if pc>rand
            u = rand;
            if (u <= 0.5)
               beta = (2*u)^(1/(GAP.mc_ec+1));
            else
               beta = (1/(2*(1-u)))^(1/(GAP.mc_ec+1));
            end
            c1_temp=0.5*((1+beta)*p1(gene)+(1-beta)*p2(gene));
            c2_temp=0.5*((1-beta)*p1(gene)+(1+beta)*p2(gene));
            c1(gene)=generepair(c1_temp,GAP);
            c2(gene)=generepair(c2_temp,GAP);
            swap=-swap;
            a=0.5*(1+swap);
            b=0.5*(1-swap);
         else
            c1(gene)= p1(gene)*a+p2(gene)*b;
            c2(gene)= p2(gene)*a+p1(gene)*b;
         end
      end
      
   end % loop over the number of chromosomes
 
% scalar simple blend based crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% fc       =  fraction of chromosomes to be crossed over
% c1       =  genes of child 1
% c2       =  genes of child 2
% GAP      =  genetic algorithm parameters
function [c1,c2]=scalar_simple_blend_crossover(p1,p2,chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               center=0.5*(p1(index)+p2(index));
               delta=2.0*(rand-0.5)*abs(p1(index)-p2(index));
               c1(index)=generepair(center+delta,GAP);
               c2(index)=generepair(center-delta,GAP);
            else
               if rand>0.5
                  c1(index)=p2(index);
                  c2(index)=p1(index);
               end
            end
         end
         
      end
              
   end % loop over the number of chromosomes
   
% vector simple blend based crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% c1       =  genes of child 1
% c2       =  genes of child 2
% GAP      =  genetic algorithm parameters
function [c1,c2]=vector_simple_blend_crossover(p1,p2,chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
         u=rand;        
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               center=0.5*(p1(index)+p2(index));
               delta=2.0*(u-0.5)*abs(p1(index)-p2(index));
               c1(index)=generepair(center+delta,GAP);
               c2(index)=generepair(center-delta,GAP);
            end
         end
         
      end
              
   end % loop over the number of chromosomes
      
   
% scalar rectangular distribution based crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% GAP      =  genetic algorithm parameters
%
% Outputs:
% c1       =  genes of child 1
% c2       =  genes of child 2
function [c1,c2]=scalar_rect_dist_crossover(p1,p2,chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;  
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               width=GAP.rdcmnw+(GAP.rdcmxw-GAP.rdcmnw)*(rand)^5;
               w1=rand>0.5;
               w2=rand>0.5;
               c1min=p1(index)*w1+(1-w1)*p2(index)-width/2.0;
               c2min=p1(index)*w2+(1-w2)*p2(index)-width/2.0;
               c1(index)=generepair(c1min+rand*width,GAP);
               c2(index)=generepair(c2min+rand*width,GAP);
            else
               if rand>0.5
                  c1(index)=p2(index);
                  c2(index)=p1(index);
               end
            end
         end
         
      end
              
   end % loop over the number of chromosomes
   
% vector rectangular distribution based crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% GAP      =  genetic algorithm parameters
%
% Outputs:
% c1       =  genes of child 1
% c2       =  genes of child 2

function [c1,c2]=vector_rect_dist_crossover(p1,p2,chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;  
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         width=GAP.rdcmnw+(GAP.rdcmxw-GAP.rdcmnw)*(rand)^5;        
        
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               c1min=p1(index)-width/2.0;
               c2max=p2(index)+width/2.0;
               c1(index)=generepair(c1min+rand*width,GAP);
               c2(index)=generepair(c2max-rand*width,GAP);
            end
         end
         
      end
              
   end % loop over the number of chromosomes

% scalar simulated binary crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% GAP      =  genetic algorithm parameters
%
% Outputs:
% c1       =  genes of child 1
% c2       =  genes of child 2

function [c1,c2]=scalar_simulated_binary_crossover(p1,p2, ...
                                                   chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;  
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
   
     
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
                % determine beta
                u = rand;
                if (u <= 0.5)
                   beta = (2*u)^(1/(GAP.mc_ec+1));
                else
                   beta = (1/(2*(1-u)))^(1/(GAP.mc_ec+1));
               end
               c1_temp=0.5*((1+beta)*p1(index)+(1-beta)*p2(index));
               c2_temp=0.5*((1-beta)*p1(index)+(1+beta)*p2(index));
               c1(index)=generepair(c1_temp,GAP);
               c2(index)=generepair(c2_temp,GAP);
            end
         end
         
      end
              
   end % loop over the number of chromosomes  
         
   
% vector simulated binary crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% GAP      =  genetic algorithm parameters
%
% Outputs:
% c1       =  genes of child 1
% c2       =  genes of child 2

function [c1,c2]=vector_simulated_binary_crossover(p1,p2,chrom_id, ...
                                                   type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;  
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
      
      % determine beta
      u = rand;
      if (u <= 0.5)
         beta = (2*u)^(1/(GAP.mc_ec+1));
      else
         beta = (1/(2*(1-u)))^(1/(GAP.mc_ec+1));
      end
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               c1_temp=0.5*((1+beta)*p1(index)+(1-beta)*p2(index));
               c2_temp=0.5*((1-beta)*p1(index)+(1+beta)*p2(index));
               c1(index)=generepair(c1_temp,GAP);
               c2(index)=generepair(c2_temp,GAP);
            end
         end
         
      end
              
   end % loop over the number of chromosomes  
   
   
% scalar quasi biological crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% fc       =  fraction of chromosomes to be crossed over
% c1       =  genes of child 1
% c2       =  genes of child 2
% GAP      =  genetic algorithm parameters
function [c1,c2]=scalar_quasi_biological_crossover(p1,p2,chrom_id, ...
                                                   type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
       
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1
               
               a=rand;
               oma=1-a;
               
               xp1=p1(index);
               xp2=p2(index);
               
               delta_xp1_min=max(-a*xp1,oma*(xp1-1));
               delta_xp1_max=min(a*(1-xp1),oma*xp1);
               delta_xp1=(delta_xp1_min+ ...
                          (delta_xp1_max-delta_xp1_min)*rand)*rand^2;
               xp1_a=  a*xp1+delta_xp1;
               xp1_b=oma*xp1-delta_xp1;
               
               delta_xp2_min=max(-a*xp2,oma*(xp2-1));
               delta_xp2_max=min(a*(1-xp2),oma*xp2);
               delta_xp2=(delta_xp2_min+ ...
                          (delta_xp2_max-delta_xp2_min)*rand)*rand^2;
               xp2_a=  a*xp2+delta_xp2;
               xp2_b=oma*xp2-delta_xp2;
               
               c1(index)=xp1_a+xp2_b;
               c2(index)=xp2_a+xp1_b;
               
            else
               if rand>0.5
                  c1(index)=p2(index);
                  c2(index)=p1(index);
               end
            end
         end
         
      end
              
   end % loop over the number of chromosomes
   
   
% vector quasi biological crossover
%
% Inputs:
% p1       =  genes of parent 1
% p2       =  genes of parent 2
% chrom_id =  mapping of the chromosome id of each gene
% type_id  =  type mapping of each chromosome
% c1       =  genes of child 1
% c2       =  genes of child 2
% GAP      =  genetic algorithm parameters
function [c1,c2]=vector_quasi_biological_crossover(p1,p2, ...
                                                   chrom_id,type_id,GAP)

   % initialize the children
   c1=p1;
   c2=p2;
   
   % do exchange on each chromosome (with probability fc) 
   for k=1:max(chrom_id) 
       
      % find list of genes on the given chromosome
      gene_list=find(chrom_id==k);
      N_genes=max(size(gene_list));
       
      % determine weather or not to do a crossover
      if  (GAP.mc_fc>rand)
        
         a=rand;
         oma=1-a;  
                  
         for i=1:N_genes
            index=gene_list(i);
            if type_id(index)>1

               xp1=p1(index);
               xp2=p2(index);
               
               delta_xp1_min=max(-a*xp1,oma*(xp1-1));
               delta_xp1_max=min(a*(1-xp1),oma*xp1);
               delta_xp1=delta_xp1_min+(delta_xp1_max-delta_xp1_min)*rand;
               xp1_a=  a*xp1+delta_xp1;
               xp1_b=oma*xp1-delta_xp1;
               
               delta_xp2_min=max(-a*xp2,oma*(xp2-1));
               delta_xp2_max=min(a*(1-xp2),oma*xp2);
               delta_xp2=delta_xp2_min+(delta_xp2_max-delta_xp2_min)*rand;
               xp2_a=  a*xp2+delta_xp2;
               xp2_b=oma*xp2-delta_xp2;
               
               c1(index)=xp1_a+xp2_b;
               c2(index)=xp2_a+xp1_b;                
                                
            end
         end
         
      end
              
   end % loop over the number of chromosomes
   
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