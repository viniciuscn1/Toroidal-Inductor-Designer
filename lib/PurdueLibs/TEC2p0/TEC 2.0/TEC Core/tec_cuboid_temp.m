function [Tmn,Tpk]=tec_cuboid_temp(C,T)

% tec_cuboid_temp calculates the steady-state spatial mean and peak
% tmperature within a cupoidal region
%
% [TAmn,TApk] = tec_cuboidal_temp(C,T)
%
% Inputs:
% C        = structure of cuboid parameters 
%            see tec_cuboid.m for documentation
% T        = vector of nodal temperatures (K)
%
% Outputs:
% TAmn     = spatial mean temperature of cuboid (K)
% TApk     = spatial peak temperature of cuboid (K)
%
% Internal:
% Q0x      = heat transfer rate into x=0 surface (W)
% Q0y      = heat transfer rate into y=0 surface (W)
% Q0z      = heat transfer rate into z=0 surface (W)
% Qlx      = heat transfer rate into x=lx surface (W)
% Qly      = heat transfer rate into y=ly surface (W)
% Qlz      = heat transfer rate into z=lz surface (W)
% c1x      = x-axis spatial temperature coefficient 1 (K/m)
% c1y      = y-axis spatial temperature coefficient 1 (K/m)
% c1z      = z-axis spatial temperature coefficient 1 (K/m)
% c2x      = x-axis spatial temperature coefficient 2(K/m^2)
% c2y      = y-axis spatial temperature coefficient 2(K/m^2)
% c2z      = z-axis spatial temperature coefficient 2(K/m^2)
% lx2      = C.lx^2 (m^2)
% ly2      = C.ly^2 (m^2)
% lz2      = C.lz^2 (m^2)
% c0       = shared spatial temperature coefficient (K)
% xe       = x-axis extremum temperature point condidate (m)
% ye       = y-axis extremum temperature point condidate (m)
% ze       = z-axis extremum temperature point condidate (m)
% Tex      = extremum temperature component in x-axis (K)
% Tey      = extremum temperature component in y-axis (K)
% Tez      = extremem temperature component in z-axis (K)
%
% Written by:
% S.D. Sudhoff                               
% Purdue University
% Electrical Engineering Building
% 465 Northwestern Avenue
% West Lafayette, IN 47907-2035
% sudhoff@purdue.edu
% 765-497-7648                               
  
% compute mean temperature
Tmn=T(C.nmn);
   
% compute heat flows    
if (C.n0x>0)
   Q0x=(T(C.n0x)-T(C.ncx))*C.Gx;
else
   if (C.n0x==0)
      Q0x=-T(C.ncx)*C.Gx;
   else
      Q0x=0;
   end
end
if (C.n0y>0)
   Q0y=(T(C.n0y)-T(C.ncy))*C.Gy;
else 
   if (C.n0y==0)
      Q0y=-T(C.ncy)*C.Gy;
   else
      Q0y=0;
   end
end
if (C.n0z>0)
   Q0z=(T(C.n0z)-T(C.ncz))*C.Gz;
else
   if (C.n0z==0) 
      Q0z=-T(C.ncz)*C.Gz;
   else
      Q0z=0;
   end
end
if (C.nlx>0)
   Qlx=(T(C.ncx)-T(C.nlx))*C.Gx;
else
   if (C.nlx==0)
      Qlx=T(C.ncx)*C.Gx;
   else
      Qlx=0;
   end
end
if (C.nly>0)
   Qly=(T(C.ncy)-T(C.nly))*C.Gy;
else
   if (C.nly==0) 
      Qly=T(C.ncy)*C.Gy;
   else
      Qly=0;
   end
end
if (C.nlz>0)
   Qlz=(T(C.ncz)-T(C.nlz))*C.Gz;
else
   if (C.nlz==0) 
      Qlz=T(C.ncz)*C.Gz;
   else
      Qlz=0.0;
   end
end

% compute spatial temperature coefficients
c1x=-Q0x./(C.kx*C.Sx);
c1y=-Q0y./(C.ky*C.Sy);
c1z=-Q0z./(C.kz*C.Sz);
  
c2x=(Q0x-Qlx)/(2*C.kx*C.V);
c2y=(Q0y-Qly)/(2*C.ky*C.V);
c2z=(Q0z-Qlz)/(2*C.kz*C.V);
   
lx2=C.lx^2;
ly2=C.ly^2;
lz2=C.lz^2;
c0=Tmn-(c2x*lx2+c2y*ly2+c2z*lz2)/3 - ...
       (c1x*C.lx+c1y*C.ly+c1z*C.lz)/2;
   
% compute extremum temperature for x-axis
Talx=c2x*lx2+c1x*C.lx;
Tex=max(0,Talx);
if (c2x~=0)
   xe=-c1x/(2*c2x);
   if (xe>=0)&&(xe<=C.lx)
      Tex=max([0 Talx -0.25*c1x^2/c2x]);
   end
end
   
% compute extremum temperature for y-axis
Taly=c2y*ly2+c1y*C.ly;
Tey=max(0,Taly);
if (c2y~=0)
   ye=-c1y/(2*c2y);
   if (ye>=0)&&(ye<=C.ly)
      Tey=max([0 Taly -0.25*c1y^2/c2y]);
   end
end

% compute extremum temperature for z-axis
Talz=c2z*lz2+c1z*C.lz;
Tez=max(0,Talz);
if (c2z~=0)
   ze=-c1z/(2*c2z);
   if (ze>=0)&&(ze<=C.lz)
      Tez=max([0 Talz -0.25*c1z^2/c2z]);
   end
end
   
% compute the peak temperatuer in the sample
Tpk=Tex+Tey+Tez+c0;

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
