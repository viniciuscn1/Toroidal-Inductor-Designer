function tec_display_material_list(TEC,vm)
% tec_display_material_list
%
% [TEC]     = tec_display_material_list(TEC)
% [TEC]     = tec_display_material_list(TEC,vm)
%
% Inputs:
% vm        = verbose mode if vm=1
% TEC       = TEC data structure (see tec_init for documentation)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

for i=1:TEC.ml.nm
   disp(['Material index = ' num2str(i) ' Description = ' TEC.ml.des{i}]);
   if (nargin==2)
      if vm==1
         disp(['   kx  = ' num2str(TEC.ml.kx(i)) ' W/m^2']);
         disp(['   ky  = ' num2str(TEC.ml.ky(i)) ' W/m^2']);
         disp(['   kz  = ' num2str(TEC.ml.kz(i)) ' W/m^2']);
         disp(['   c   = ' num2str(TEC.ml.c(i)) ' J/(kg*K)']);
         disp(['   row = ' num2str(TEC.ml.row(i)) ' kg/m^3']);
      end
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
