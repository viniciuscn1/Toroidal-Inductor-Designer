function [penalty] = divcon(P,GAP)
% DIVCON  Computes a penalty function in order to ensure genetic diversity.
%
% [penalty] = divcon(P,GAP)
%
% Inputs:
% P        = population structure
% GAP      = genetic algorithm parameter structure
%
% Outputs:
%
% penalty  = output population (after diversity penalty has been
%            calculated)
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

% return if diversity control is not used
if ~GAP.dc_act
   if max(size(P))>1
      penalty=ones(1,size(P,2));    
   else
      penalty=ones(1,P.size);
   end
   return
end

% initialize the penalty
penalty=ones(1,P.size);

% for each subset population
Nregions=max(P.region);
for region=1:Nregions,
    
    sub_population_list=find(P.region==region);
    
    if GAP.dc_spc == 1   % diversity control on parameter space
        subpop = P.normgene(:,sub_population_list);
    else                 % diversity control on the fitness function space
        subpop = P.mfit(:,sub_population_list);
        subpopmx = max(subpop,[],2);
        subpopmn = min(subpop,[],2);
        subpop = (subpop-subpopmn*ones(1,size(subpop,2)))./ ...
                 ((subpopmx-subpopmn)*ones(1,size(subpop,2)));
    end
    
    switch GAP.dc_alg
        case {1}
           penalty(sub_population_list)=divcona(subpop,GAP.dc_mnt, ...
                                                GAP.dc_mxt); 
        case {2}
           penalty(sub_population_list)=divconb(subpop,GAP.dc_ntr, ...
                                                GAP.dc_mnb,GAP.dc_mxb);         
        case {3}
           penalty(sub_population_list)=divconc(subpop,GAP.dc_dc);
        case {4}
           penalty(sub_population_list)=divcond(subpop,GAP.dc_dc, ...
                                                GAP.dc_nt); 
        otherwise
           error('Requested diversity control algorithm nonexistant');
    end
   
end 

function p=divcona(x,mnt,mxt)

   % determine the number of elements
   n=size(x,2);
   
   % initialize the penalty
   p=ones(1,n);
   
   % determine the threshold
   thr=mnt+rand*(mxt-mnt);
      
   % compute the distance between every element and every other element
   d=zeros(n,n);
   for i=1:n,
      xi=x(:,i);
      for j=i+1:n,
         xj=x(:,j);
         d(i,j)=norm(xj-xi);
         d(j,i)=d(i,j);
      end
   end
   
   % compute the mean distance between points times the threshold
   mdthr=median(median(d))*thr;
   
   % determine the penalty
   for i=1:n,
       p(i)=1.0/(length(find(d(:,i)<mdthr)));
   end
         

function p=divconb(x,ntr,mnb,mxb)

   % determine the number of elements in each vector 
   % and the number of vectors
   [ne,nv]=size(x);
   
   % compute limits on the number of bins
   nbmin=nv*mnb;
   nbmax=nv*mxb;
   
   % initialize the penalty
   p=zeros(1,nv);
   
   % compute vector penalty for each element
   for e=1:ntr,
       
      % compute the number of bins
      nb=round(nbmin+(nbmax-nbmin)*rand);
   
      % compute the edges of the bins
      be=linspace(0,1+eps,nb+1);
      
      % now use a hashing function
      hf=randperm(ne);
      xh=rem(hf*x,1);
      
      % get the count and bin index for each bin
      [c,bi]=histc(xh,be);
      
      % update penalty for each vector
      for v=1:nv,
          penalty=1/c(bi(v));        % penalty based on count
          if penalty>p(v)
             p(v)=penalty;
          end
      end

   end % loop over each element
   
function p=divconc(x,dc)
      
   % determine the number of vectors and number of elements
   [~,nv]=size(x);
   
   % compute the distance between every element and every other element
   d=zeros(nv,nv);
   for i=1:nv,
      xi=x(:,i);
      for j=i+1:nv,
         xj=x(:,j);
         d(i,j)=norm(xj-xi,inf);
         d(j,i)=d(i,j);
      end
   end

   % determine the count
   c(nv)=0.0;
   for i=1:nv,
       c(i)=sum(exp(-d(i,:)/dc));
   end
   
   % compute the penalty
   p=1./c;
     
function p=divcond(x,dc,nt)
      
   % determine the number of vectors and number of elements
   [~,nv]=size(x);
   
   % initialize the penalty
   p=ones(1,nv);
   
   % compute the distance between every vector and a subset of vector
   % the subset never includes the vector itself
   d=zeros(nv,nt-1);
   for i=1:nv,
      xi=x(:,i);
      k=1;
      while k<nt
         j=ceil(nv*rand);
         if j==0
            j=1; 
         end
         if j~=i
            xj=x(:,j);
            d(i,k)=norm(xj-xi,inf);
            k=k+1;
         end
      end
   end
  
   % determine the count
   c(nv)=0.0;
   for i=1:nv,
       c(i)=1+sum(exp(-d(i,:)/dc))*nv/nt;
   end
 
   % compute the penalty
   p=1./c;
   
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
