function [CP]=conductor_catalog(n)
% conductor_catalog assigns the appropriate conductor parameters
%                   to a structure of parameters.
%
% Call:
% CP=conductor_catalog(n)
%
% Inputs:
% n          = conductor identification number
%              n=1 - copper
%              n=2 - aluminum
%              n=3 - silver
%              n=4 - gold
%
% Outputs:
% CP         = structure of conductor parameters
%  CP.desc   = conductor description
%  CP.t0     = temperature at which conductivity characterized in C
%  CP.sigma0 = conductivity at t0 (1/(Ohm-m))
%  CP.sigmac = same as CP.sigma0 (backwards compatibility)
%  CP.alpha  = temperature coefficient of resisitivity (1/K)
%              assuming the relationship 
%              [(1/sigma)| at t] = [(1/sigmac)| at t0]*(1+alpha(t-t0));
%  CP.jmax   = recommended maximum current density (A/m^2) 
%  CP.row    = mass density (kg/m^3)
%  CP.k      = thermal conductivity (W/(K-m))
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
%
% References:
% Fink, D. G. (ed.), and H. W. Beaty (ed.), Standard Handbook 
% for Electrical Engineers. 13th ed. New York: McGraw-Hill, Inc., 1993.
%
% Lide, D. E. (ed.), CRC Handbook of Chemistry and Physics.
% 84th ed. Boca Raton, FL: CRC Press, 2003.
%
% Serway, R. A., Principles of Physics. 2nd ed. Fort Worth,
% TX:  Saunders College Pub., 1998.
%
% Stauffer, H. B., Engineer's Guide to the National Electric
% Code. Boston, MA: Jones and Bartlett Publishers, Inc., 2008.
%
% Young, H. D., University Physics.  7th ed. Reading, MA: 
% Addison-Wesley, 1992.
%
% Note: This code was taken form the Power Magnetic Material Toolbox by
%       S.D. Sudhoff

switch n
   case 1
      CP.desc='Copper';
      CP.sigma0=5.959e7; % [Lide, sec. 12]
      CP.t0=20;
      CP.alpha=3.93e-3;  % [Fink, pg. 4-9]
      CP.jmax=7.6025e6;  % calculated using 10 AWG at 90 C (194 F) which 
                         % has a max ampacity of 40 A [Stauffer, Table 3.3]
      CP.row=8890;       % @ 20 C (pure copper, rolled, forged, or drawn 
                         % and then annealed) [Fink, pg. 4-3]
      CP.k=385;          % [Young, Table 15-5] 
   case 2
      CP.desc='Aluminum';
      CP.sigma0=3.774e7; % [Lide, sec. 12]
      CP.t0=20;
      CP.alpha=3.94e-3;  % [Fink, pg. 4-10]
      CP.jmax=6.6522e6;  % calculated using 10 AWG at 90 C (194 F) which 
                         % has a max ampacity of 35 A [Stauffer, Table 3.3]
      CP.row=2705;       % @ 20 C (commercially hard-drawn) [Fink, pg. 4-3]
      CP.k=205;          % [Young, Table 15-5] 
   case 3
      CP.desc='Silver';
      CP.sigma0=6.301e7; % [Lide, sec. 12]
      CP.t0=20;
      CP.alpha=3.8e-3;   % [Serway, pg. 602]
      CP.jmax=NaN;
      CP.row=10499;      % [Fink, Table 4-34]
      CP.k=406;          % [Young, Table 15-5] 
   case 4
      CP.desc='Gold';
      CP.sigma0=4.517e7; % [Lide, sec. 12]
      CP.t0=20;
      CP.alpha=3.4e-3;   % [Serway, pg. 602]
      CP.jmax=NaN;
      CP.row=19301;      % [Fink, Table 4-34]
      CP.k=314;          % [Young, Table 15-5]      
   otherwise
      CP.desc='Bad';
      CP.sigma0=NaN;
      CP.t0=NaN;
      CP.alpha=NaN;
      CP.jmax=NaN;
      CP.row=NaN;
      CP.k=NaN;
              
end

   CP.sigmac=CP.sigma0;
          
end

%  Copyright 2013 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of the Power Magnetic Material Toolbox.
% 
%  The Power Magnetic Material Toolbox is free software: 
%  you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  The Power Magnetic Material Toolbox is distributed 
%  in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with the Power Magnetic Material Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.