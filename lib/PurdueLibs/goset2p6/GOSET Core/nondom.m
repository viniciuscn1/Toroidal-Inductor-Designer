function [nd] = nondom(f,t)
% NONDOM  Finds the set of non-dominated solutions
%
% [nd] = nondom(f,t)
%
% Inputs:
% f        = a (number of objectives) by (number of solutions) 
%            matrix of objective function values
% t        = type flag.  1 indicates larger value of objective values are 
%            better; 0 indicates smaller values are better.
%
% Outputs:
% nd       = a row vector with dimension equal to number of solutions whose
%            elements are 1 if the solutions are nondominated, or 0 if they
%            are dominated
%
% Approach 3 from K. Deb, Multi-objective Optimization Using Evolutionary 
%    Algorithms, John Wiley and Sons, Inc., 2001.   pp. 38-39

f = (2*t - 1)*f;

[f1,ind1] = sort(f(1,:));   % account for type and sort according
                            % to first objective function
f1 = [f1; f(2:end,ind1)];   % append the rest of the obj function values
f1 = fliplr(f1);            % flip so they are in descending order
ind1 = fliplr(ind1);

% take care of case where f1 values are identical
[~,c] = size(f1);
for k = 1:(c-1)
    
    if (f1(1,k) == f1(1,k+1))
        i1 = find(f1(2:end,k) > f1(2:end,k+1));
        if isempty(i1)
            f1(:,[k (k+1)]) = f1(:,[(k+1) k]);
            ind1([k (k+1)]) = ind1([(k+1) k]);
        end
    end
    
end

nd_ind = front(f1);         % call recursive function first time
nd = zeros(1,size(f,2));    % calculate output
nd(ind1(nd_ind)) = ones(1,length(nd_ind));



function nd_ind = front(P)

% implement recursive function 'front()' used in this approach

[~,I] = size(P);

if (I == 1)
    
    nd_ind = 1;
    
else
    
    ind_T = 1:floor(I/2);
    ind_B = (floor(I/2)+1):I;
    
    ndT = front(P(:,ind_T));
    ndB = front(P(:,ind_B));
    
    addtoT = [];
    for i = 1:length(ndB);
        
        isdominated1 = all(P(:,ind_B(ndB(i)))*ones(1,length(ndT)) <= ...
                       P(:,ind_T(ndT)));
        isdominated2 = any(P(:,ind_B(ndB(i)))*ones(1,length(ndT)) < ...
                       P(:,ind_T(ndT)));
        isdominated = any(isdominated1.*isdominated2);
        
        if ~isdominated
            addtoT = [addtoT ind_B(ndB(i))];
        end
        
    end
    
    nd_ind = [ind_T(ndT) addtoT];
    
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

        