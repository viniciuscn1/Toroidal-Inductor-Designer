function [geneval,fit] = despick(bI,f)

% DESPICK  A function to pick and obtain fitness and gene values of a
%          particular non-dominated solution from the Pareto Optimal
%          front. It also plots the Pareto-Optimal front with the chosen
%          solution highlighted. This works only with 2 objective
%          optimization.
%
%
% [geneval,fit] = despick(bI,f)
% [geneval]     = despick(bI,f)
%
%
% Inputs:
% bI       = Best individuals (For the two objective case,
%            each column is a non-dominated individual(sorted according to 
%            the specified objective).
% f        = Best fitness/fitnesses(Each column is a set of fitness values
%            of the individuals corresponding to bI)
%
% Note:  bI and f can be obtained as outputs from GAOPTIMIZE
%
% Outputs:
% gene_val = The gene parameter values of the chosen non-dominated design
% fit      = The fitness value/values of the chosen non-dominated design

[r1,c1] = size(f);
if r1 ~= 2
    error('This code works only for 2 objective optimization')
end
a1 = int2str(c1);
desnum = ...
input(['Please enter non-dominated solution number less than ',a1,' :']) ;
if desnum > c1
    error('Wrong input')
end

obj1_srt = f(1,:);
obj2_srt = f(2,:);
figure
plot(obj1_srt,obj2_srt,'-xb','Markersize',6);
hold on
plot(obj1_srt(desnum),obj2_srt(desnum), ...
     'pr','MarkerSize',10,'MarkerFaceColor','r');
% Labeling
xlabel('Objective 1');
ylabel('Objective 2');
title('Pareto Optimal Front') ;
legend('All solutions',['Design ',int2str(desnum)]) ;

switch nargout
    case 2
        geneval = bI(:,desnum);
        fit  = f(:,desnum);
    case 1
        geneval = bI(:,desnum);
        fit = [];
    case 0
        geneval = [];
        fit = [];
    otherwise
        error('Wrong number of output arguments')
end
end

%  Script by Harish Suryanarayana
%  Modified Sep 7th, 2010

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
%  License along with GOSET(x.x). If not, see 
%  <http://www.gnu.org/licenses/>.