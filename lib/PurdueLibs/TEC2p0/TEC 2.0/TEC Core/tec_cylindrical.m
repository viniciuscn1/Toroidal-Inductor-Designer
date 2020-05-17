function [TEC,CE] = tec_cylindrical(TEC,mc,ri,ro,lz,theta, ...
                                    nir,n0z,ncr,ncz,nmn, ...
                                    nor,nlz,p,gp)                                                
% tec_cylindrical  Adds a hollow cylindrical element 
%                  to a thermal equivalent circuit
%
% [TEC]      = tec_cylindrical(TEC,mc,ri,ro,lz,theta,nir,n0z,ncr,ncz, ...
%                              nmn,nor,nlz,p,gp)
%
% Inputs: 
% TEC      = TEC data structure (see tec_init for documentation)
% mc       = material code
%            Note: thermal conuctivity in x direction is used as radial
%                  direction
% ri       = inner radius of cylinder (m)
% ro       = outer radius of cylinder (m)
% lz       = length in z-direction (m)
% theta    = angle spanned by cyllinder (m)
% nir      = node at r=ri (use -1 for open circuit)
% n0z      = node number at z=0 (use -1 for open circuit)
% ncr      = node number at center of radial axis of 'T' circuit 
% ncz      = node number at center of z-axis 'T' circuit
% nmn      = node number for mean temperature of cuboid 
% nor      = node at r=ro (use -1 for open circuit)
% nlz      = node number at z=lz (use -1 for open circuit)
% p        = nominal power into cuboid (W)
% gp       = slope of power dissipation w.r.t. mean temperature (W/K)
%
% Outputs:
% CE       = structure of cyclindrical element parameters
%  CE.V    = volume (m^3)
%  CE.M    = mass (kg)
%  CE.C    = heat capacity (J/K)
%  CE.Gx   = thermal conductance from edge to center in x-direction (W/K)
%  CE.Gy   = thermal conductance from edge to center in y-direction (W/K)
%  CE.Gz   = thermal conductance from edge to center in z-direction (W/K)
%  CE.Sx   = area perpendicular to x-axis (m^2)
%  CE.Sy   = area perpendicular to y-axis (m^2)
%  CE.Sz   = area perpendicular to z-axis (m^2)
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

% compute volume, mass, thermal capacitance,dimensions
CE.Sri=theta*ri*lz;
CE.Sro=theta*ro*lz;
CE.Sz=0.5*theta*(ro^2-ri^2);
CE.V=CE.Sz*lz;
CE.M=CE.V*TEC.ml.row(mc);
CE.C=CE.M*TEC.ml.c(mc);
CE.ri=ri;
CE.ro=ro;
CE.lz=lz;
CE.theta=theta;

% record node numbers
CE.nir=nir;
CE.ncr=ncr;
CE.nor=nor;
CE.n0z=n0z;
CE.ncz=ncz;
CE.nlz=nlz;
CE.nmn=nmn;

% check kr=kx=ky
if TEC.ml.kx(mc)~=TEC.ml.ky(mc)
    error('kr = kx = ky, but kx = %.3f and ky = %.3f',...
        TEC.ml.kx(mc),TEC.ml.ky(mc))
end

% record thermal conductivities
CE.kr = TEC.ml.kx(mc);
CE.kz = TEC.ml.kz(mc);

% compute thermal conductances
ro2     = ro^2;
ri2     = ri^2;
ro2mri2 = ro2-ri2;
ro2pri2 = ro2+ri2;
lnrori  = log(ro/ri);
ro2ri2  = ro2*ri2;
twokrthetal = 2.0*CE.kr*theta*lz;
CE.Gir  = twokrthetal/(2*ro2*lnrori/ro2mri2-1);
CE.Gor  = twokrthetal/(1-2*ri2*lnrori/ro2mri2);
CE.Gtr  = -2*twokrthetal*ro2mri2/ ...
        (ro2pri2-4*ro2ri2*lnrori/ro2mri2);
CE.Gz   = 2*CE.kz*CE.Sz/lz;

% put in branches
if (nir>=0)TEC=tec_g_branch(TEC,nir,ncr,CE.Gir); end;
if (nor>=0)TEC=tec_g_branch(TEC,nor,ncr,CE.Gor); end;
if (n0z>=0)TEC=tec_g_branch(TEC,n0z,ncz,CE.Gz); end;
if (nlz>=0)TEC=tec_g_branch(TEC,nlz,ncz,CE.Gz); end;
TEC=tec_g_branch(TEC,ncr,nmn,CE.Gtr);
TEC=tec_g_branch(TEC,ncz,nmn,(-3*CE.Gz));
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
