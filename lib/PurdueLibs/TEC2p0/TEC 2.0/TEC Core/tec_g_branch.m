function [TEC] = tec_g_branch(TEC,pn,nn,gb)
% tec_g_branch  Adds a conductance branch to a thermal equivalent circuit
%
% [TEC]     = tec_g_branch(TEC,pn,nn,gb)
%
% Inputs:
% TEC      = TEC data structure {see tec_init for details}
% pn       = postive node
% nn       = negative node
% gb       = thermal conductance of branch (W/K)
%
% Outputs:
% TEC      = Updated TEC data structure
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

if pn>0
   TEC.Gn(pn,pn)=TEC.Gn(pn,pn)+gb;
end
if nn>0
   TEC.Gn(nn,nn)=TEC.Gn(nn,nn)+gb;
end
if (pn>0)&&(nn>0)
   TEC.Gn(pn,nn)=TEC.Gn(pn,nn)-gb;
   TEC.Gn(nn,pn)=TEC.Gn(nn,pn)-gb;
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
