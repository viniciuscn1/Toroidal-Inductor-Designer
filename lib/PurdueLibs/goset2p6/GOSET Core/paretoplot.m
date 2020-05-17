function paretoplot(P,GAP,region)
% PARETOPLOT Plots two or three objective function fitness 
%            against each other
%
% paretoplot(P,GAP,[region])
%
% Inputs:
% P        = chromosome population
% GAP      = genetic algorithm parameter structure
% region   = an optional integer argument which, if included
%            will only plot the chromosomes in a given region
%
% Written by:
% S.D. Sudhoff
% Purdue University
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% check to make sure all chromosomes have been evaluated
if min(P.eval)<1
   error('ERROR in distplot: population not fully evaluated');
end

% find number of objetives to plot
if ~isempty(~GAP.pp_list)

   % number of objectives to plot 
   nop=size(GAP.pp_list,2);                 

   % make a list of population list
   if nargin==2
      pop_list=1:P.size;
   else
      if nargin==3
         pop_list=find(P.region==region);
     else
         error('Error in pareto2dplot: invalid number of arguments'); 
      end
   end
   
   % 2-d pareto optimal plot  
   if (nop==2)
   
      % determine the objectives we will be plotting
      obj1=GAP.pp_list(1);
      obj2=GAP.pp_list(2);
         
      % check the settings
      if ~((GAP.pp_style(obj1)==0)||(GAP.pp_style(obj1)==1))
         error('Unknown Pareto Plot Style');
      end
      if ~((GAP.pp_style(obj2)==0)||(GAP.pp_style(obj2)==1))
         error('Unknown Pareto Plot Style');
      end
      if ~((GAP.pp_sign(obj1)==1)||(GAP.pp_sign(obj1)==-1))
         error('Pareto Plot Sign Parameter Must Be +1 or -1');
      end
      if ~((GAP.pp_sign(obj2)==1)||(GAP.pp_sign(obj2)==-1))
         error('Pareto Plot Sign Parameter Must Be +1 or -1');
      end
         
      % find the conditioned fitness
      f1=GAP.pp_sign(obj1)*P.mfit(obj1,pop_list);
      f2=GAP.pp_sign(obj2)*P.mfit(obj2,pop_list);
      if GAP.pp_style(obj1)==0
         f1=GAP.pp_sign(obj1)*log10(f1);
      end
      if GAP.pp_style(obj2)==0
         f2=GAP.pp_sign(obj2)*log10(f2);
      end
         
      % determine dominate and non-dominated individuals
      nd=nondom(P.mfit,1);
      ndi=find(nd==1);
      di=find(nd==0);
        
      % plot individuals
      plot(f1(di),f2(di),'bx');
      hold on;
      plot(f1(ndi),f2(ndi),'ro');
      hold off;
      grid;
      xlabel(GAP.pp_xl);
      ylabel(GAP.pp_yl);
      title(GAP.pp_title);
      if ~isempty(GAP.pp_axis)
         axis(GAP.pp_axis);
      end
      drawnow;
     
   end % nop=2 
   
   % 3-d pareto optimal plot  
   if (nop==3)
   
      % determine the objectives we will be plotting
      obj1=GAP.pp_list(1);
      obj2=GAP.pp_list(2);
      obj3=GAP.pp_list(3);
      
      % check the settings
      if ~((GAP.pp_style(obj1)==0)||(GAP.pp_style(obj1)==1))
         error('Unknown Pareto Plot Style');
      end
      if ~((GAP.pp_style(obj2)==0)||(GAP.pp_style(obj2)==1))
         error('Unknown Pareto Plot Style');
      end
      if ~((GAP.pp_style(obj3)==0)||(GAP.pp_style(obj3)==1))
         error('Unknown Pareto Plot Style');
      end
      if ~((GAP.pp_sign(obj1)==1)||(GAP.pp_sign(obj1)==-1))
         error('Pareto Plot Sign Parameter Must Be +1 or -1');
      end
      if ~((GAP.pp_sign(obj2)==1)||(GAP.pp_sign(obj2)==-1))
         error('Pareto Plot Sign Parameter Must Be +1 or -1');
      end
      if ~((GAP.pp_sign(obj3)==1)||(GAP.pp_sign(obj3)==-1))
         error('Pareto Plot Sign Parameter Must Be +1 or -1');
      end
         
      % find the conditioned fitness
      f1=GAP.pp_sign(obj1)*P.mfit(obj1,pop_list);
      f2=GAP.pp_sign(obj2)*P.mfit(obj2,pop_list);
      f3=GAP.pp_sign(obj3)*P.mfit(obj3,pop_list);
      if GAP.pp_style(obj1)==0
         f1=GAP.pp_sign(obj1)*log10(f1);
      end
      if GAP.pp_style(obj2)==0
         f2=GAP.pp_sign(obj2)*log10(f2);
      end
      if GAP.pp_style(obj3)==0
         f3=GAP.pp_sign(obj3)*log10(f3);
      end
         
      % determine dominate and non-dominated individuals
      nd=nondom(P.mfit,1);
      ndi=find(nd==1);
      di=find(nd==0);
        
      set(gcf,'Userdata',view);  % Store current view angle
      
      % plot individuals
      plot3(f1(di),f2(di),f3(di),'bx');
      hold on;
      plot3(f1(ndi),f2(ndi),f3(ndi),'ro');
      hold off;
      grid;
      xlabel(GAP.pp_xl);
      ylabel(GAP.pp_yl);
      zlabel(GAP.pp_zl);
      title(GAP.pp_title);
      if ~isempty(GAP.pp_axis)
         axis(GAP.pp_axis);
      end
      
      view(get(gcf,'Userdata'));    % Restore 3D view angle
      
      drawnow;
     
   end % nop=3    
   
   
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