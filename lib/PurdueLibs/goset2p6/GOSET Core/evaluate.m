function [mfit,es,une] = evaluate(P,GAP,cne,D)
% EVALUATE  Evaluates the fitness of a population of chromosomes.
%           Note that this routine only updates those memebers of
%           the population whos fitness has not already been evaluated
%           Observe that all genes of a populaton member are passed, as
%           well as the region number for that member of the population
%
% [mfit,es,une] = evaluate(P,GAP,cne,[D])
%
% Inputs:
% P              = current population of chromosomes (a structure)
% GAP            = a structure with the genetic algorith parameters
% cne            = current number of evaluations performed thus far
% D              = an optional argument or structure with all data needed  
%                  to evaluate the function
%
% Outputs:
% mfit           = multiobjective fitness (number of objectives by size of 
%                  population)
% es             = evaluation status of each member of population
% une            = updated number of evaluations
%
% Internal:
% extra_data_present = true if extra data present
% mfit           = fitness values
% update_indices = population members needing evaluation updates
% nupdates       = number of updates to make
% fithandle      = fitness handle
% nepg0          = smallest number of evaluations in a parallel group
% nepg           = number of evaluations in each group
% se             = array of indices to start evaluation
% ee             = array if indices to end evaluation
% temp           = cell array with fitness values for each group
% gene           = cell array with gene values for each group
% age            = cell array with age values for each group
% omfit          = cell array of old mfit values for each group
% region         = cell array with region values for each group
% grp_ind        = cell array with index values for each group
%
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

% check on the number of input argument
if ~((nargin==3)||(nargin==4))
   error('Error in Evaluate: Wrong number of input arguments');
end

% determine if the extra data needs to be passed to the evaluation routine 
extra_data_present=(nargin==4);
if extra_data_present
   extra_data_present=size(D,1)>0;
end

% initialize the fitness matrix
mfit=P.mfit;

% find the elements whos indices need to be updated
if GAP.ev_are
   update_indices=1:P.size;
else
   update_indices=find(P.eval==0);
end

% find the number of updates to make  
nupdates=length(update_indices);

% see if there are any updates to make
if nupdates==0
   une=cne;
   es=ones(1,size(mfit,2));
   return;
end

% single evaluation w/o parallel computing---------------------------------
if (~GAP.ev_bev)&&(~GAP.ev_pp)
  
   if extra_data_present
      if GAP.ev_ssd 
         for i=1:nupdates
            j=update_indices(i);
            mfit(:,j)=feval(P.fithandle,P.gene(:,j),P.age(j), ...
                            P.mfit(:,j),P.region(j),D);
         end
      else
         for i=1:nupdates
            j=update_indices(i);
            mfit(:,j)=feval(P.fithandle,P.gene(:,j),D);
         end
      end
   else
      if GAP.ev_ssd
         for i=1:nupdates
             j=update_indices(i);
             mfit(:,j)=feval(P.fithandle,P.gene(:,j),P.age(j), ...
                             P.mfit(:,j),P.region(j));
         end
      else
         for i=1:nupdates
            j=update_indices(i);
            mfit(:,j)=feval(P.fithandle,P.gene(:,j));
         end
      end
   end
    
   
   % update evaluation status (since everything is now evaluated)
   es=ones(1,size(mfit,2));

   % statistics
   une=cne+nupdates;
   
   return
    
end %----------------------------------------------------------------------    
    
    

% block evaluation w/o parallel computing----------------------------------
if GAP.ev_bev&&(~GAP.ev_pp)
   
    if extra_data_present
      if GAP.ev_ssd 
         mfit(:,update_indices)=feval(P.fithandle, ...
                                      P.gene(:,update_indices), ...
                                      P.age(update_indices), ...
                                      P.mfit(:,update_indices), ...
                                      P.region(:,update_indices),D);
      else
         mfit(:,update_indices)=feval(P.fithandle, ...
                                      P.gene(:,update_indices),D);
      end
   else
      if GAP.ev_ssd
         mfit(:,update_indices)=feval(P.fithandle, ...
                                      P.gene(:,update_indices), ...
                                      P.age(update_indices), ...
                                      P.mfit(:,update_indices), ...
                                      P.region(:,update_indices));
      else
         mfit(:,update_indices)=feval(P.fithandle, ...
                                      P.gene(:,update_indices));
      end	
   end    
   
   % update evaluation status (since everything is now evaluated)
   es=ones(1,size(mfit,2));

   % statistics
   une=cne+nupdates;
   
   return
      
end %----------------------------------------------------------------------
     

% single evaluation with parallel computing--------------------------------
if (~GAP.ev_bev)&&GAP.ev_pp   

   
   % divide up population data into parrallel groups
   fithandle=P.fithandle;
   nepg0=floor(nupdates/GAP.ev_npg);
   nepg=nepg0*ones(GAP.ev_npg,1);
   nepg(1:nupdates-GAP.ev_npg*nepg0)=nepg(1:nupdates-GAP.ev_npg*nepg0)+1;
   ee=cumsum(nepg);
   se=ee-nepg+1;
   temp=cell(GAP.ev_npg,1);
   gene=cell(GAP.ev_npg,1);
   age=cell(GAP.ev_npg,1);
   omfit=cell(GAP.ev_npg,1);
   region=cell(GAP.ev_npg,1);
   grp_ind=cell(GAP.ev_npg,1);
   for i=1:GAP.ev_npg
       grp_ind{i}=update_indices(se(i):ee(i));
       gene{i}=P.gene(:,grp_ind{i});
       age{i}=P.age(grp_ind{i});
       omfit{i}=P.mfit(:,grp_ind{i});
       region{i}=P.region(grp_ind{i});
   end
         
   % single evaluation
   if extra_data_present
      if GAP.ev_ssd 
  
         parfor i=1:GAP.ev_npg
            for j=1:nepg(i)
               temp{i}(:,j)=feval(fithandle,gene{i}(:,j),age{i}(j), ...
                                  omfit{i}(:,j),region{i}(j),D); 
            end
         end
                  
      else
          
         parfor i=1:GAP.ev_npg
            for j=1:nepg(i)
               temp{i}(:,j)=feval(fithandle,gene{i}(:,j),D); 
            end
         end
                  
      end
   else
      if GAP.ev_ssd
        
         parfor i=1:GAP.ev_npg
            for j=1:nepg(i)
               temp{i}(:,j)=feval(fithandle,gene{i}(:,j),age{i}(j), ...
                                  omfit{i}(:,j),region{i}(j)); 
            end
         end
         
      
      else

         parfor i=1:GAP.ev_npg
            for j=1:nepg(i)
               temp{i}(:,j)=feval(fithandle,gene{i}(:,j)); 
            end
         end
      
      end	
   end    

   % populate mfit
   for i=1:GAP.ev_npg
      mfit(:,grp_ind{i})=temp{i}; 
   end
      
   % update evaluation status (since everything is now evaluated)
   es=ones(1,size(mfit,2));

   % statistics
   une=cne+nupdates;
   
   return
      
end %----------------------------------------------------------------------
      
% block evaluation with parallel computing---------------------------------
if (GAP.ev_bev)&&GAP.ev_pp      
         
    % divide up population data into parrallel groups
   fithandle=P.fithandle;
   nepg0=floor(nupdates/GAP.ev_npg);
   nepg=nepg0*ones(GAP.ev_npg,1);
   nepg(1:nupdates-GAP.ev_npg*nepg0)=nepg(1:nupdates-GAP.ev_npg*nepg0)+1;
   ee=cumsum(nepg);
   se=ee-nepg+1;
   temp=cell(GAP.ev_npg,1);
   gene=cell(GAP.ev_npg,1);
   age=cell(GAP.ev_npg,1);
   omfit=cell(GAP.ev_npg,1);
   region=cell(GAP.ev_npg,1);
   grp_ind=cell(GAP.ev_npg,1);
   for i=1:GAP.ev_npg
       grp_ind{i}=update_indices(se(i):ee(i));
       gene{i}=P.gene(:,grp_ind{i});
       age{i}=P.age(grp_ind{i});
       omfit{i}=P.mfit(:,grp_ind{i});
       region{i}=P.region(grp_ind{i});
   end
         
   % single evaluation
   if extra_data_present
      if GAP.ev_ssd 
  
         parfor i=1:GAP.ev_npg
            temp{i}=feval(fithandle,gene{i},age{i}(j), ...
                          omfit{i}(:,j),region{i}(j),D); 
         end
                          
      else
          
         parfor i=1:GAP.ev_npg
            temp{i}=feval(fithandle,gene{i},D); 
         end
                         
      end
   else
      if GAP.ev_ssd
        
         parfor i=1:GAP.ev_npg
             temp{i}=feval(fithandle,gene{i},age{i}(j), ...
                           omfit{i}(:,j),region{i}(j)); 
         end
              
      else
          
         parfor i=1:GAP.ev_npg
            temp{i}=feval(fithandle,gene{i}); 
         end
      
      end	
   end    

   % populate mfit
   for i=1:GAP.ev_npg
      mfit(:,grp_ind{i})=temp{i}; 
   end
        
   % update evaluation status (since everything is now evaluated)
   es=ones(1,size(mfit,2));

   % statistics
   une=cne+nupdates;

   return
      
end %----------------------------------------------------------------------

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