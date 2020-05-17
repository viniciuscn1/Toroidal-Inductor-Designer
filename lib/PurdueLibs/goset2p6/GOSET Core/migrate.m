function [O] = migrate(N,GAP,cg)
% MIGRATE performs a migration operation
%
% [O] = migrate(N,GAP,cg)
%
% Inputs:
% N        = population before migration
% GAP      = structure containing genetic algorithm parameters
%            (see gap_default)
% cg       = current generation number
%
% Outputs:
% O        = population after migration (a structure)
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


% iniitialize migrated population
O=N;

% determine whether to use migration
if (GAP.mg_nreg>1)&&(GAP.mg_pmig>0)

   % determine number of regions
   nregions=max(N.region);
    
   % initialize the number of generation since last migration
   persistent next_migration
   if isempty(next_migration) 
      next_migration=cg+ceil(0.5+(0.5+rand)*GAP.mg_tmig);
   elseif cg == 2
      next_migration = 2;
   end

   % update generation and perform migration when 
   % appropriate
   if (next_migration==cg)
      % compute number of migrations 
      N_migrations=O.size*GAP.mg_pmig;
   
      for i=1:N_migrations,
       
         % pick an individual to migrate 
         individual=ceil(max([rand eps])*N.size); 
      
         % pick out a new region
         new_region=ceil(max([rand eps])*nregions);
      
         % assign the new region
         O.region(individual)=new_region;
         O.eval(individual)=0;
      
      end
      
      % reset the generation number
      next_migration=next_migration + ceil(0.5+(0.5+rand)*GAP.mg_tmig);   
  end
      
end % migration being used

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
