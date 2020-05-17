function reportplot(GAP,GAS,Pk)
% REPORTPLOT  Plots the distribution of the genes of the chromosomes
%             and the fitness history 
%
% reportplot(GAP,GAS,Pk)
%
% Inputs:
% GAP      = Genetic Algorithm Parameters
%            (see gap_default for definitions)
% GAS      = A structure containting the statistics of the GA
%            (see GAOPTIMIZE.M for definitions)
% Pk       = chromosome population

if (GAS.cg == 1)||(GAS.cg == GAP.fp_ngen)||(rem(GAS.cg,GAP.rp_gbr)< 1.0e-8)

if (GAP.rp_lvl == 1)
    
   for i=1:length(GAP.op_list),
    
      % find the objective 
      obj = GAP.op_list(i); 
      
      % name of plot
      objplotname = ['Gene Distribution / Objective ',num2str(obj), ...
                     ' Fitness'];
       
      figHandles = allchild(0);
      figH = findobj(figHandles, 'flat', 'Name', objplotname,'Tag', ...
                     'reportplot');
      if ~isempty(figH)
         figure(figH)
      else
         figH = figure('NumberTitle', 'off', 'Name', objplotname,...
                       'Tag', 'reportplot');
      end

      subplot(2,1,1);
      distplot(figH,Pk,obj,GAP);
      subplot(2,1,2);
      
      if ~((GAP.op_style(obj)==1)||(GAP.op_style(obj)==0))
         error('Invalid Value of Plotting Parameter GAP.op_style');
      end
      if ~((GAP.op_sign(obj)==1)||(GAP.op_sign(obj)==-1))
         error('Invalid Value of Plotting Parameter GAP.op_sign');
      end
      
      % compute modified plotting values
      bestfit   = GAP.op_sign(obj)*GAS.bestfit(obj,1:GAS.cg);
   
      meanfit   = GAP.op_sign(obj)*GAS.meanfit(obj,1:GAS.cg);   
      medianfit = GAP.op_sign(obj)*GAS.medianfit(obj,1:GAS.cg);
      if GAP.op_style(obj) == 0
            bestfit   = GAP.op_sign(obj)*log10(bestfit);
            meanfit   = GAP.op_sign(obj)*log10(meanfit);   
            medianfit = GAP.op_sign(obj)*log10(medianfit);
      end
      
      % do the plot
      plot(1:GAS.cg,medianfit,'g', ...
           1:GAS.cg,meanfit,'r', ...
           1:GAS.cg,bestfit,'b');
      set(gca,'XLIM',[1 GAP.fp_ngen]);
      if GAP.op_style(obj)==1
         if GAP.op_sign(obj)==1 
            ylabel('f:B(b) Md(g) Mn(r)');
         else
            ylabel('-f:B(b) Md(g) Mn(r)');
         end
      else
         if GAP.op_sign==1
            ylabel('log10(f):B(b) Md(g) Mn(r)');
         else
            ylabel('-log10(-f):B(b) Md(g) Mn(r)');
         end
      end
      xlabel('Generation');
    
   end
   
   if ~isempty(GAP.pp_list)
      figHandles = allchild(0);
      figH=findobj(figHandles, 'flat', 'Name', 'Pareto Plot','Tag', ...
                   'reportplot');
      if ~isempty(figH)
         figure(figH)         
      else
         figH = figure('NumberTitle', 'off', 'Name', 'Pareto Plot',...
                       'Tag', 'reportplot'); 
         % Set the default view angle
         view(size(GAP.pp_list,2));  
      end
      paretoplot(Pk,GAP);
      drawnow;
   end
   
   % custom report plot
   if ~isempty(GAP.rp_crh)
       feval(GAP.rp_crh,Pk,GAS);
   end
    
end

if (GAP.rp_lvl~=-1)
    disp(strcat(['Statistics for generation ' num2str(GAS.cg)]));
    disp(strcat(['Best fitness = ' num2str(GAS.bestfit(:,GAS.cg)')]));
    disp(strcat(['Mean fitness = ' num2str(GAS.meanfit(:,GAS.cg)')]));
    disp(strcat(['Median fitness = ' num2str(GAS.medianfit(:,GAS.cg)')]));
    disp(strcat(['Number of evaluations = ' num2str(GAS.ne(GAS.cg))]));
    disp(' ')
end

end

drawnow;    % Refresh the plot

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