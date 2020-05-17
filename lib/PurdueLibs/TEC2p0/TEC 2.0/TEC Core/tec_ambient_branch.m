function [TEC] = tec_ambient_branch(TEC,n,S,h,Ta)

% tec_ambient  Adds a connention to a node to an ambient temperature source
%
% [TEC]     = tec_ambient(TEC,n,S,h,Ta)
%
% Inputs:
% TEC      = TEC data structure (see tec_init for documentation)
% n        = node number
% S        = contact area (m^2)
% h        = heat transfer coefficient (W/(K*m^2))
% Ta       = ambient temperature source (K)
%
% Outputs:
% TEC      = Updated TEC data structure
%
% Internal:
% g        = conductance to ambient source (W/K)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

g=h*S;
TEC.Gn(n,n)=TEC.Gn(n,n)+g;
TEC.Pn(n)=g*Ta;

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



