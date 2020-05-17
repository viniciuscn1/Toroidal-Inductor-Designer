function [Tmn,Tpk]=tec_onedim_temp(O,T)
%
% tec_onedim_temp  Find the mean and peak temperature within a 
%                  one-dimension heat flow element 
%
% Inputs:
% O        = structure of one-dimensional heat flow element parameters     
% T        = vector of nodal temperatures (K)
%
% Outputs:
% TAmn     = spatial mean temperature of element (K)
% TApk     = spatial peak temperature of elememt (K)
% 
% Internal:
% T0x      = temperature at x=0 node (K)
% Tlx      = temperature at x=lx node (K)
% Tcx      = temperature at center node (K)
% Tmn      = temperature at mean node (K)
% P        = power dissipation within element (W)
% p        = power dissipation density within element (W/m^3)
% xex      = extemum point of termperature (m)
% Tex      = extrumum temperature (K)
%                   
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu  


Tcx=T(O.ncx);
Tmn=T(O.nmn);

if O.n0x>0
   T0x=T(O.n0x);
else
   if (O.n0x==0)
      T0x=0.0;
   else
      T0x=Tcx;
   end
end
if O.nlx>0
   Tlx=T(O.nlx);
else
   if (O.nlx==0)
      Tlx=0.0;
   else
      Tlx=Tcx;
   end
end

P=(Tcx-Tmn)*3*O.Gx;
p=P/(O.Sx*O.lx);

Tpk=max(T0x,Tlx);
if (p~=0)
   xex=((Tlx-T0x)/O.lx+p*O.lx/(2*O.kx))*O.kx/p;
   if (xex>=0)&&(xex<=O.lx)
      Tex=(0.5*O.kx/p)*((Tlx-T0x)/O.lx+0.5*p*O.lx/O.kx)^2+T0x;
      Tpk=max([T0x Tlx Tex]);
   end
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
