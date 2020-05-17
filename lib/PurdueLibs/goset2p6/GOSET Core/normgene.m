function [P]=normgene(P)
% NORMGENE Updates the normalized genes based on raw genes.
%          Only the population members who have not been evaluated
%          are updated.
%
% P=normgene(P)
%
% Inputs:
% P = population structure
%
% Outputs:
% P= population structure with updated normalized genes
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
plist = find(P.eval==0);

% normalize the integer variables to the specified range
for type=1:3,
   tlist=find(P.type==type);
   ntlist=length(tlist);
   switch type
       case {1,2}
          for i=1:ntlist,
             gindex=tlist(i); 
             a=P.min(gindex);
             b=P.max(gindex)-a;
             if b>0
                P.normgene(gindex,plist)=(P.gene(gindex,plist)-a)/b;
             else
                P.normgene(gindex,plist)=0;
             end
          end
       case {3}
          for i=1:ntlist,
             gindex=tlist(i); 
             a=log(P.min(gindex));
             b=log(P.max(gindex))-a;
             if b>0
                P.normgene(gindex,plist)=(log(P.gene(gindex,plist))-a)/b;
             else
                P.normgene(gindex,plist)=0;
             end
          end
   end
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
