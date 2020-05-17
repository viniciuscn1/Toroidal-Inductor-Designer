function [rgene]=generepair(gene,GAP)
% GENEREPAIR Repairs illegal gene values.
%
% rgene=generepair(gene)
%
% Inputs:
% gene = an individual, vector, or array of normalized gene values
%
% Outputs:
% rgene = repaired gene values
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

% find list of individuals to evaluate
rgene=gene;
repair_list1=find(gene>1);
repair_list2=find(gene<0);

switch GAP.gr_alg
    case 1 % hard limit
       rgene(repair_list1)=1;
       rgene(repair_list2)=0;
    case 2 % ring mapping
       rgene(repair_list1)=mod(rgene(repair_list1),1);
       rgene(repair_list2)=mod(rgene(repair_list2),1);
    otherwise
       error('Unknown gene repair algorithm');
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
