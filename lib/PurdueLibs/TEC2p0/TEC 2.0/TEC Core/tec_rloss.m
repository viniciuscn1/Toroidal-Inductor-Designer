function [P] = tec_rloss(C,T,J,Vc)

% tec_rloss  Computes resistive loss in an element
%
% [TEC]    = tec_cuboid(TEC,mx,lx,ly,lz,n0x,n0y,n0z,nlx,nly,nlz,nmn,p)
%
% Inputs:
% TEC      = TEC data structure (see tec_init for documentation)
% C        = Structure of conductor paramaters
%  C.sigma = Conductivity at temperature of C.T0 (S)
%  C.alpha = Temperature coefficient of conductivity
%  C.T0    = Temperature at which nominal conductivity specified (K)
% J        = current density (A/m^2)
% Vc       = volume of actual conductor within element (m^3)
%
% Outputs:
% P        = resistive power dissipation in element (W)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

sigma=C.sigma0/(1+C.alpha*(T-C.T0));
P=Vc*J^2/sigma;

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


