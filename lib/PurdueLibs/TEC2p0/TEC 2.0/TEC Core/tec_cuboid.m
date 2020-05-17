function [TEC,CB] = tec_cuboid(TEC,mc,lx,ly,lz, ...
                               n0x,n0y,n0z,ncx,ncy,ncz,nmn, ...
                               nlx,nly,nlz,p,gp)
% tec_cuboid  Adds a cuboid to a thermal equivalent circuit
%
% [TEC]     = tec_cuboid(TEC,mx,lx,ly,lz,n0x,n0y,n0z, ...
%                        ncx,ncy,ncz,nmn,nlx,nly,nlz,p,gp)
%
% Inputs:
% TEC      = TEC data structure (see tec_init for documentation)
% mc       = material code 
% lx       = length in x-direction (m)
% ly       = length in y-direction (m)
% lz       = length in z-direction (m)
% n0x      = node number at x=0 (use -1 for open circuit)
% n0y      = node number at y=0 (use -1 for open circuit)
% n0z      = node number at z=0 (use -1 for open circuit)
% ncx      = node number at center of x-axis 'T' circuit
% ncy      = node number at center of y-axis 'T' circuit
% ncz      = node number at center of z-axis 'T' circuit
% nmn      = node number for mean temperature of cuboid 
% nlx      = node number at x=lx (use -1 for open circuit)
% nly      = node number at y=ly (use -1 for open circuit)
% nlz      = node number at z=lz (use -1 for open circuit)
% p        = nominal power into cuboid (W)
% gp       = slope of power dissipation w.r.t. mean temperature (W/K)
%
% Outputs:
% CB       = structure of cupoid parameters
%  CB.V    = volume (m^3)
%  CB.M    = mass (kg)
%  CB.C    = heat capacity (J/K)
%  CB.Gx   = thermal conductance from edge to center in x-direction (W/K)
%  CB.Gy   = thermal conductance from edge to center in y-direction (W/K)
%  CB.Gz   = thermal conductance from edge to center in z-direction (W/K)
%  CB.Sx   = area perpendicular to x-axis (m^2)
%  CB.Sy   = area perpendicular to y-axis (m^2)
%  CB.Sz   = area perpendicular to z-axis (m^2)
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
CB.V=lx*ly*lz;
CB.M=CB.V*TEC.ml.row(mc);
CB.C=CB.M*TEC.ml.c(mc);
CB.Sx=ly*lz;
CB.Sy=lx*lz;
CB.Sz=lx*ly;
CB.lx=lx;
CB.ly=ly;
CB.lz=lz;

% record node numbers
CB.n0x=n0x;
CB.ncx=ncx;
CB.nlx=nlx;
CB.n0y=n0y;
CB.ncy=ncy;
CB.nly=nly;
CB.n0z=n0z;
CB.ncz=ncz;
CB.nlz=nlz;
CB.nmn=nmn;

% record thermal conductivities
CB.kx=TEC.ml.kx(mc);
CB.ky=TEC.ml.ky(mc);
CB.kz=TEC.ml.kz(mc);

% compute thermal conductances
CB.Gx=2*CB.kx*CB.Sx/lx;
CB.Gy=2*CB.ky*CB.Sy/ly;
CB.Gz=2*CB.kz*CB.Sz/lz;

% put in branches
if (n0x>=0)TEC=tec_g_branch(TEC,n0x,ncx,CB.Gx); end;
if (nlx>=0)TEC=tec_g_branch(TEC,nlx,ncx,CB.Gx); end;
if (n0y>=0)TEC=tec_g_branch(TEC,n0y,ncy,CB.Gy); end;
if (nly>=0)TEC=tec_g_branch(TEC,nly,ncy,CB.Gy); end;
if (n0z>=0)TEC=tec_g_branch(TEC,n0z,ncz,CB.Gz); end;
if (nlz>=0)TEC=tec_g_branch(TEC,nlz,ncz,CB.Gz); end;
TEC=tec_g_branch(TEC,ncx,nmn,(-3*CB.Gx));
TEC=tec_g_branch(TEC,ncy,nmn,(-3*CB.Gy));
TEC=tec_g_branch(TEC,ncz,nmn,(-3*CB.Gz));
TEC=tec_std_branch(TEC,nmn,0 ,0 ,gp,nmn,0 ,p);

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

