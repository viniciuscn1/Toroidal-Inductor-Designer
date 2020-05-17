function [TEC,OD] = tec_onedim(TEC,mc,lx,Sx, ...
                               n0x,ncx,nmn,nlx,p,gp)
%
% tec_onedim  Adds a one-dimensional heat flow element to a 
%             thermal equivalent circuit
%
% [TEC]     = tec_onedim(TEC,mc,lx,Sx,n0x,ncx,)
%
% Inputs:
% TEC      = TEC data structure (see tec_init for documentation)
% mc       = material code 
% lx       = length in x-direction (m)
% Sx       = cross section perpindicular to x-direction (m)
% n0x      = node number at x=0 (use -1 for open circuit)
% ncx      = node number at center of x-axis 'T' circuit
% nmn      = node number for mean temperature of cuboid 
% nlx      = node number at x=lx (use -1 for open circuit)
% p        = nominal power into cuboid (W)
% gp       = slope of power dissipation w.r.t. mean temperature (W/K)
%
% Outputs:
% OD       = structure of one-dimensional heat flow element parameters
%  OD.V    = volume (m^3)
%  OD.M    = mass (kg)
%  OD.C    = heat capacity (J/K)
%  OD.Gx   = thermal conductance from edge to center in x-direction (W/K)
%  OD.Sx   = area perpendicular to x-axis (m^2)
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

% compute volume, mass, thermal capacitance
OD.V=lx*Sx;
OD.M=OD.V*TEC.ml.row(mc);
OD.C=OD.M*TEC.ml.c(mc);
OD.Sx=Sx;
OD.lx=lx;

% record node numbers
OD.n0x=n0x;
OD.ncx=ncx;
OD.nlx=nlx;
OD.nmn=nmn;

% record thermal conductivities
OD.kx=TEC.ml.kx(mc);

% compute thermal conductances
OD.Gx=2*OD.kx*OD.Sx/lx;

% put in branches
if (n0x>=0)TEC=tec_g_branch(TEC,n0x,ncx,OD.Gx); end;
if (nlx>=0)TEC=tec_g_branch(TEC,nlx,ncx,OD.Gx); end;
TEC=tec_g_branch(TEC,ncx,nmn,(-3*OD.Gx));
TEC=tec_std_branch(TEC,nmn,0,0,gp,nmn,0,p);

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

