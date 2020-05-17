function [TEC] = tec_init(nn,emnb)
% tec_init   Initializes an Thermal Equivalent Circuit analysis.
%            All fields are initialized to zero
%
% [TEC]     = tec_init(nn,emnb)
% [TEC]     = tec_init(nn)
%
% Inputs:
% nn        = number of nodes
% emnb      = estimate of maximum number of branches
%
% Outputs:
% TEC          = initialized TEC data sturcture
%  TEC.nn      = number of non-zero nodes.  Node 0 is ground.
%  TEC.Gn      = network thermal conductance matrix (TEC.nn by TEC.nn) (W/K)
%  TEC.Pn      = network power vector (TEC.nn by 1) (W)
%  Material List Stucture
%  TEC.ml      = material list structure
%   TEC.ml.nm  = number of materials in material library
%   The following vectors are TEC.ml.nm by 1
%   TEC.ml.des = material description
%   TEC.ml.k   = thermal conductivity (K/W)
%   TEC.ml.c   = specific heat capacity (J / (K*kg))
%   TEC.ml.row = mass density (kg / m^3)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initialize nodal quantitities
TEC.nn=nn;
TEC.Gn=zeros(nn,nn);
TEC.Pn=zeros(nn,1);

% initialize materials
TEC.ml.nm=0;
[TEC,~]=tec_material(TEC,'generic copper',385,390,8890);
[TEC,~]=tec_material(TEC,'generic aluminum',205,910,2705);
[TEC,~]=tec_material(TEC,'generic M19',16.7,469,7402); 
[TEC,~]=tec_material(TEC,'generic M36',18.8,465,7018); 
[TEC,~]=tec_material(TEC,'generic M43',20.9,461,7291);
[TEC,~]=tec_material(TEC,'generic M47',37.7,446,7585);
[TEC,~]=tec_material(TEC,'generic Hiperco 50',29.8,418,7845);
[TEC,~]=tec_material(TEC,'air at 300K',0.0263,1007,1.161);   
[TEC,~]=tec_material(TEC,'generic magnet wire insulation',0.4,1500,1400);

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
