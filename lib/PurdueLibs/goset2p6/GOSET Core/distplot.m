function distplot(fignum,P,objective,GAP,region)
% DISTPLOT Plots the distribution of the genes of the chromosomes
%
% distplot(fignum,P,objective,GAP,[region])
%
% Inputs:
% fignum   = figure number
% P        = chromosome population
% objective= objective number to show in plot
% GAP      = genetic algorithm parameter vector
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
%
% Modificatiosn
% removed drawnow command.  Scott Sudhoff 3/26/2010

% check to make sure all chromosomes have been evaluated
for i=1:P.size
   if P.eval(i)~=1
       error('ERROR in distplot: population not fully evaluated');
   end
end

% make a list of population list
if nargin==4
   pop_list=1:P.size;
   pop_size=P.size;
else
   if nargin==5
      pop_list=find(P.region==region);
      pop_size=max(size(pop_list));
   else
      error('Error in distplot: invalid number of arguments'); 
   end
end

if GAP.dp_type==1

   % select members of the population for plotting
   pop_list=pop_list(randperm(pop_size));
   pop_list=pop_list(1:min(pop_size,GAP.dp_np)); 
   pop_size=max(size(pop_list));

   % sort the chromosomes by fitediness
   [genefit,geneind]=sort(P.mfit(objective,pop_list),2,'ascend');

   % plot chromosomes
   figure(fignum);
   plot([0 0 P.ngenes+0.5 P.ngenes+0.5 0],[0 1 1 0 0],'k-');
   hold on;
   fmin=min(P.mfit(objective,:));
   fmax=max(P.mfit(objective,:));
   fmean=mean(P.mfit(objective,:));
   fmnmx=abs(fmin-fmax);
   fmdmx=abs(fmean-fmax);
   
   for i=1:pop_size
       rank=i;
       ind =pop_list(geneind(i));
       dist=rank/pop_size;
       colors=linecolor(rank,pop_size);
       plot((1:P.ngenes)-0.5+dist,P.normgene(:,ind)','.','Color',colors);
   end
     
   for i=1:P.ngenes,
       plot([i-0.5 i-0.5 i+0.5 i+0.5],[0 1 1 0],'k-');
   end

   axis([0.5 P.ngenes+0.5 0 1]);
   xlabels=linspace(1,P.ngenes,P.ngenes);
   h=gca;
   set(h,'XTick',xlabels);
   set(h,'XTickLabel',xlabels);
   xlabel('Parameter Number');
   ylabel('Normalized Value');
   hold off;

end

if GAP.dp_type==2
   
   % plot boundaries for the distribution of genes
   figure(fignum);
   plot([1 1 P.ngenes+1 P.ngenes+1 0],[0 1 1 0 0],'k-');
   hold on;
   for i=1:P.ngenes,
       plot([i i],[0 1],'k-');
   end

   % set up the bins
   res=GAP.dp_res;
   edges=linspace(0,1.00001,res+1);
   xvals=zeros(4*res);
   yvals=zeros(4*res);

   % draw the histograms
   for gene=1:P.ngenes
      n=histc(P.normgene(gene,pop_list),edges)/pop_size; 
      n(res)=n(res)+n(res+1);
      for bin=1:res
         x1=gene;
         x2=x1+n(bin);
         y1=(bin-1)/res;
         y2=bin/res;
         xvals(4*(bin-1)+1:4*bin)=[x1 x2 x2 x1];
         yvals(4*(bin-1)+1:4*bin)=[y1 y1 y2 y2];
      end
      fill(xvals,yvals,'b'); 
   end
   
  
   % define colors
   c=['b- ';'g- ';'r- ';'c- ';'m- ';'b-o';'g-o';'r-o';'c-o';'m-o'; ...
      'b-x';'g-x';'r-x';'c-x';'m-x';'b-*';'g-*';'r-*';'c-*';'m-*'];
   
   % show location of best genes
   x2vals=zeros(2*P.ngenes);
   y2vals=zeros(2*P.ngenes);
   if nargin==4
      for i=1:max(P.region)
         region_list=find(P.region==i);
         [rbestfit,rbestsubind]=max(P.mfit(objective,region_list));
         rbestindex=region_list(rbestsubind);
         for j=1:P.ngenes
             plot([j j+1], ...
                  [1 1]*P.normgene(j,rbestindex),c(mod(i-1,20)+1,:));
         end
      end
   else
      region_list=find(P.region==region);
      [rbestfit,rbestsubind]=max(P.mfit(objective,region_list));
      rbestindex=region_list(rbestsubind);
      for j=1:P.ngenes
          plot([j j+1],[1 1]*P.normgene(j,rbestindex),'g-');
      end  
   end
   
   % label the axis
   axis([1 P.ngenes+1 0 1]);
   xlabels=linspace(1,P.ngenes,P.ngenes);
   h=gca;
   set(h,'XTick',xlabels+0.5);
   set(h,'XTickLabel',xlabels);
   xlabel('Parameter Number');
   ylabel('Normalized Value');
   hold off;
   
end
end

function [colors]=linecolor(index,maxindex)

   f=(index-1)/(maxindex-1);
   blue = max(0,1-2*f);
   red  = max(0,2*f-1);
   green=1-red-blue;
   colors=[red green blue];
   if (min(colors)<0)||(max(colors)>1)
       keyboard
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
