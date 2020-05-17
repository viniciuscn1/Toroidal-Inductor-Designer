function [Tmn,Tpk]=tec_cylindrical_temp(C,T)

% tec_cylindrical_temp calculates the steady-state spatial mean and peak
% temperature within a cylindrical region
%
% [TAmn,TApk] = tec_cylindrical_temp(C,T)
%
% Inputs:
% C        = structure of cylindrical element parameters 
%            see tec_cylindrical.m for documentation
% T        = vector of nodal temperatures (K)
%
% Outputs:
% TAmn     = spatial mean temperature of cylindrical element (K)
% TApk     = spatial peak temperature of cylindrical element (K)
%
% Internal:
% Qir      = heat transfer rate into r=ri surface (W)
% Q0z      = heat transfer rate into z=0 surface (W)
% Qor      = heat transfer rate out of r=ro surface (W)
% Qlz      = heat transfer rate out of z=lz surface (W)
% c2r      = radial spatial temperature coefficient (r^2 term) (K/m^2)
% clr      = radial spatial temperature coefficeint (ln term) (K)
% c1z      = z-axis spatial temperature coefficient 1 (K/m)
% c2z      = z-axis spatial temperature coefficient 2(K/m^2)
% c0       = spatial temperature coefficient 0 (K)
% ro2      = intermediate variable (m^2)
% ri2      = intermediate variable (m^2)
% lz2      = intermediate variable (m^2)
% re       = radial extremum temperature point condidate (m)
% ze       = z-axis extremum temperature point condidate (m)
% Tri      = intermediate variable (K)
% Tro      = intermediate variable (K)
% Ter      = extremum temperature component in radial direction (K)
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
if (C.nir>0)
   Qir=(T(C.nir)-T(C.ncr))*C.Gir;
else
   if (C.nir==0)    
      Qir=-T(C.ncr)*C.Gir;
   else
      Qir=0;
   end
end
if (C.n0z>0)
   Q0z=(T(C.n0z)-T(C.ncz))*C.Gz;
else
   if (C.noz==0)
      Q0z=-T(C.ncz)*C.Gz;
   else
      Q0z=0;
   end
end
if (C.nor>0)
   Qor=(T(C.ncr)-T(C.nor))*C.Gor;
else
   if (C.nor==0)
      Qor=T(C.ncr)*C.Gor;
   else
      Qor=0;
   end
end
if (C.nlz>0)
   Qlz=(T(C.ncz)-T(C.nlz))*C.Gz;
else
   if (C.nlz==0) 
      Qlz=T(C.ncz)*C.Gz;
   else
      Qlz=0;
   end
end

% compute spatial distribution terms
ro2=C.ro^2;
ri2=C.ri^2;
lz2=C.lz^2;
clr=(Qor*ri2-Qir*ro2)/(2*C.kr*C.V);
c1z=-Q0z./(C.kz*C.Sz);
c2r=(Qir-Qor)/(4*C.kr*C.V);
c2z=(Q0z-Qlz)/(2*C.kz*C.V);
c0=Tmn-(0.5*c2r*(ro2+ri2)+ ...
        clr*((ro2*log(C.ro)-ri2*log(C.ri))/(ro2-ri2)-0.5)+...
        c2z*lz2/3+c1z*lz2/2);
   
% compute extremum temperature for radial term
Tri=c2r*ri2+clr*log(C.ri);
Tro=c2r*ro2+clr*log(C.ro);
Ter=max(Tri,Tro);
if (clr*c2r<0)
   re=sqrt(-0.5*clr/c2r);
   if (re>C.ri)&&(re<C.ro)
      Ter=max([Tri Tro c2r*re^2+clr*log(re)]);
   end
end
   
% compute extremum temperature for z-axis
Talz=c2z*lz2+c1z*C.lz;
Tez=max(0,Talz);
if (c2z~=0)
   ze=-c1z/(2*c2z);
   if (ze>0)&&(ze<C.lz)
      Tez=max([0 Talz -0.25*c1z^2/c2z]);
   end
end
   
% compute the peak temperatuer in the sample
Tpk=Ter+Tez+c0;

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
