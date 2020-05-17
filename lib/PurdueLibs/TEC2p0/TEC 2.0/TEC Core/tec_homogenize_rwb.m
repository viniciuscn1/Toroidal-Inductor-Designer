function [TEC,mih]=tec_homogenize_rwb(TEC,w,d,N,rc,ti,mic,mii,mia,des)

%  tec_homogenize_rwb computes a homogenized represenation of a w (x-axis) 
%                     by d (y-axis) winding bundle. Determines effective 
%                     thermal conductivities in x-,y-, and z-(length) axis
%                     as well as effective density and effective specific 
%                     heat capacity
%
% [TEC,mih]=tec_homogenize_rwb(TEC,w,d,N,rc,ti,mic,mii,mia,des)
%
% Inputs:
% TEC      = thermal equivalent circuit structure   
% w        = width of conductor bundle (m) 
% d        = depth of conductor bundle (m)
% N        = number of conductors
% rc       = conductor radius (m)
% ti       = insulation thickness (m)
% mic      = material index of conductor 
% mii      = material index of insulator
% mia      = material index of air/potting material
% des      = description of homogenized material
%
% Outputs:
% TEC      = updated thermal equivalent circuit structure
%
% Internal:
% ac       = total condutor cross section (m^2)
% ai       = total insulator cross section (m^2)
% aa       = total air/potting cross section (m^2)
% zeta     = winding aspect ratio
% dc       = depth of conductor region (y-axis) (m)
% di       = depth of insulator region 'bot L' (y-axis) (m)
% da       = depth of air/potting region 'bot L'(y-axis)(m)
% kc       = thermal conductivity of the conductor (W/(m*K))
% ki       = thermal conductivity of insulation (W/(m*K))
% ka       = thermal conductivity of air/potting (W/(m*K))
% rowc     = mass density of conductor (W/(m*K))
% rowi     = mass density of insulator (W/(m*K))
% rowa     = mass density of air (W/(m*K))
% cc       = specific heat capacity of conductor (J/(K*kg))
% ci       = specific heat capacity of insulator (J/(K*kg))
% ca       = specific heat cap. of air/potting (J/(K*kg))
% kxy      = effective thermal conductivity in xy 
%            directions (W/(m*K))
% kz       = effective thermal conductivity in z
%            direction (W/(m*K))
% rowe     = effective density (kg/m^3)
% ce       = effective specific heat capacity
%
% Written by
% Scott D. Sudhoff                               
% Purdue University
% Electrical Engineering Building
% 465 Northwestern Avenue
% West Lafayette, IN 47907-2035
% sudhoff@ecn.purdue.edu
% 765-497-7648

% look up material properties
kc=TEC.ml.kx(mic);
ki=TEC.ml.kx(mii);
ka=TEC.ml.kx(mia);
cc=TEC.ml.c(mic);
ci=TEC.ml.c(mii);
ca=TEC.ml.c(mia);
rowc=TEC.ml.row(mic);
rowi=TEC.ml.row(mii);
rowa=TEC.ml.row(mia);

% compute material areas
ac=pi*N.*rc.^2;
ai=pi*N.*((rc+ti).^2-rc.^2);
aa=w.*d-ac-ai;

% get aspect ratio
zeta=w./d;

% compute dimensions of regions
dc = sqrt(ac./zeta);
di = sqrt(dc.^2+ai./zeta)-dc;
da = (d-dc-di);

% compute thermal conductivities
kxye=1./(1./kc+di./(ki.*dc)+da./(ka.*dc)) + ...
     1./((dc+di)./(ki.*di)+da./(ka.*di)) + ...
     1./(d./(ka.*da));
kze=(ac.*kc+ai.*ki+aa.*ka)./(w.*d);
rowe=(ac.*rowc+ai.*rowi+aa.*rowa)./(w.*d);
ce=(ac.*rowc.*cc+ai.*rowi.*ci+aa.*rowa.*ca)./(w.*d.*rowe);

% add to material list
mih=TEC.ml.nm+1;
TEC.ml.nm=mih;
TEC.ml.kx(mih)=kxye;
TEC.ml.ky(mih)=kxye;
TEC.ml.kz(mih)=kze;
TEC.ml.c(mih)=ce;
TEC.ml.row(mih)=rowe;
TEC.ml.des{mih}=des;

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

  
