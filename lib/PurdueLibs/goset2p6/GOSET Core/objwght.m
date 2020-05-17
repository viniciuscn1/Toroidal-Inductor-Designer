function owv=objwght(GAP)
% OBJWGHT  Creates an objective weight vector for use in multi-objective
%          optimization.
%
% owv=objwght(GAP)
%
% Inputs:
% GAP = genetic algorithm parameter structure
%
% Outputs:
% owv = the objective weight vector
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

% single objective case
if (GAP.fp_nobj==1)
   owv=1;
   return;
end

% multi-objective case (but concentrating on a single objective)
if (GAP.fp_obj>0)
   owv=zeros(1,GAP.fp_nobj);
   owv(1,GAP.fp_obj)=1;
   return;
end

% multi-objective case (with multiobjective optimization)
if (GAP.fp_obj==0)
   owv=abs(randn(1,GAP.fp_nobj));
   owv=owv/(sum(owv)+eps);
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
