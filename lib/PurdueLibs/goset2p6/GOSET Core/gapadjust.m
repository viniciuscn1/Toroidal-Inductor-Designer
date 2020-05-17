function [GAP]= gapadjust(GAPic,GAS)
% GAPADJUST Sets parameters of genetic algorithm as a function of
%           generation number
%
% [GAP] = gapadjust(GAS,GAPic)
%
% Inputs:
% GAS        = Structure of Genetic Algorithm Statistics
% GAPic      = Initial Conditions for Genetic Algorithm Parameters
%
% Outputs:
% GAP        = Genetic Algorithm Parameters
%              Open gapdefault.m with a text editor for a full 
%              description of all parameters and default values of
%              all parameters in the GAP structure
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

GAP=GAPic;

beta=GAS.cg/GAP.fp_ngen;
alpha=1-beta;

GAP.mt_ptgm=GAPic.mt_ptgm0*alpha+GAPic.mt_ptgmf*beta;
GAP.mt_prgm=GAPic.mt_prgm0*alpha+GAPic.mt_prgmf*beta;
GAP.mt_srgm=GAPic.mt_srgm0*alpha+GAPic.mt_srgmf*beta;
GAP.mt_pagm=GAPic.mt_pagm0*alpha+GAPic.mt_pagmf*beta;
GAP.mt_sagm=GAPic.mt_sagm0*alpha+GAPic.mt_sagmf*beta;
GAP.mt_prvm=GAPic.mt_prvm0*alpha+GAPic.mt_prvmf*beta;
GAP.mt_srvm=GAPic.mt_srvm0*alpha+GAPic.mt_srvmf*beta;
GAP.mt_pavm=GAPic.mt_pavm0*alpha+GAPic.mt_pavmf*beta;
GAP.mt_savm=GAPic.mt_savm0*alpha+GAPic.mt_savmf*beta;
GAP.mt_pigm=GAPic.mt_pigm0*alpha+GAPic.mt_savmf*beta;
GAP.el_fpe = GAPic.el_fpe0*alpha+ GAPic.el_fpef*beta;

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
