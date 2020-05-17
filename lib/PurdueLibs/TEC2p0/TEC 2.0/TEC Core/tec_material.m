function [TEC,mi] = tec_material(TEC,des,kx,ky,kz,c,row)
% tec_material  Adds a material to a thermal equivalent ciruit
%
% [TEC,mi]  = tec_material(TEC,des,k,c,row)
% [TEC,mi]  = tec_material(TEC,des,kx,ky,kz,c,row)
%
% Inputs:
% TEC      = TEC data structure (see tec_init for documentation)
% des      = description.  A string.
% k        = thermal conductivity (W/(m*K))
% c        = specific heat capacity (J/(kg*K))
% row      = density (kg/m^3)
%
% Outputs:
% TEC      = Updated TEC data structure
% mi       = material index.  
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% reorder inputs
if nargin==5
   row=kz;
   c=ky;
   ky=kx;
   kz=kx;
end

% assign fields
mi=TEC.ml.nm+1;
TEC.ml.nm=mi;
TEC.ml.des{mi}=des;
TEC.ml.kx(mi)=kx;
TEC.ml.ky(mi)=ky;
TEC.ml.kz(mi)=kz;
TEC.ml.c(mi)=c;
TEC.ml.row(mi)=row;

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
