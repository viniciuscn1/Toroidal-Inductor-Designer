function [GAP]= gapdefault(nobj,obj,npop,ngen)

% GAPDEFAULT assigns default values to the genetic algorithm parameters
%             used in GAP.
%
% [GAP] = gapdefault
% [GAP] = gapdefault(nobj)
% [GAP] = gapdefault(nobj,obj,npop,ngen)
%
% Inputs:
% nobj       = Number of objectives (default is 1)
% obj        = Objective to optimize (set obj=0 for multiobjective 
%              optimization, which is default if nobj>1)
% npop       = Nominal population size (default is 100) 
% ngen       = Number of generations (default is 100)
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

% Set the number of objective if not specified
switch nargin
    case 0
       GAP.fp_nobj=1;                 % number of objectives
       GAP.fp_obj=1;                  % objective to optimize
       GAP.fp_ngen=100;               % number of genenerations to evolve 
       GAP.fp_ipop=100;               % initial population size   
       GAP.fp_npop=100;               % normal population size
    case 1
       GAP.fp_nobj=nobj;              % number of objectives
       if nobj>1
          GAP.fp_obj=0;               % objective to optimize
       else
          GAP.fp_obj=1;
       end
       GAP.fp_ngen=100;               % number of genenerations to evolve 
       GAP.fp_ipop=100;               % initial population size
       GAP.fp_npop=100;               % normal population size
    case 4
       GAP.fp_nobj=nobj;              % number of objectives
       GAP.fp_obj=obj;                % objective to optimize
       GAP.fp_ngen=ngen;              % number of genenerations to evolve 
       GAP.fp_ipop=npop;              % initial population size    
       GAP.fp_npop=npop;              % normal population size
    otherwise,
       error('Inproper call to gapdefault');       
end
                             
% Diversity control parameters
GAP.dc_act=1;                  % activate diversity control 
                               % (1=active; 0=non-active)
GAP.dc_alg=4;                  % diversity control algorithm
GAP.dc_spc=1;                  % diversity control space 
                               %    1 = Paramater (or gene) space 
                               %    2 = Fitness function space
GAP.dc_mnt=0.02;               % minimum threshold (for algorithm 1)
GAP.dc_mxt=0.1;                % maximum threshold (for algorithm 1)
GAP.dc_ntr=3;                  % number of trials (for algorithm 2)
GAP.dc_mnb=0.5;                % minumum number of bins relative to 
                               % population size (for algorithm 2)
GAP.dc_mxb=2.0;                % maximum number of bins relative to 
                               % population size (for algorithm 2)
GAP.dc_dc=1.0e-3;              % diversity control distance constraint
                               % (algorithms 3 and 4)
GAP.dc_nt=50;                  % diversity control test population size 
                               % (algorithm 4)                               

% Scaling parameters
GAP.sc_alg=1;     % algorithm
                  %    0 = none
                  %    1 = offset so minimum fitness is zero
                  %    2 = linear scaling (most fit individual GAP.sc_kln 
                  %        more likely to be selected than average fit)
                  %    3 = linear scaling (most fit individual GAP.sc_kln  
                  %        more likely to be selected than median fit)
                  %    4 = linear scaling (most fit individual GAP.sc_kln 
                  %        more likely to be selected than least fit)
                  %    5 = sigma truncation
                  %    6 = quadtratic scaling
GAP.sc_kln=10.0;  % scaling factor for linear scaling algorithms
GAP.sc_cst=2.0;   % scaling constant for sigma truncation
GAP.sc_kmxq=10.0; % scaling factor for quadratic scaling (most fit ind.
                  %    GAP.sc_kmxq more likely selected than median fit)
GAP.sc_kmnq=0.01; % scaling factor for quadratic scaling (most fit ind.
                  %    GAP.sc_kmnq more likely selected than median fit)

% Selection algorithm parameters
GAP.sl_alg=2;                  % algorithm
                               %    1 = roulette wheel
                               %    2 = tournament
                               %    3 = custom
GAP.sl_nts=4;                  % number of individuals in a tournament
GAP.sl_cah=[];                 % custom algorithm handle
                               % parent_list=f(region,size,age,mfit,fit)
                               % where
                               %    parent_list=indices of those 
                               %                individuals to be parents
                               %    region=region number
                               %    size=number of individuals to be in 
                               %         parent list
                               %    age=vector describing ages of 
                               %        individuals in region
                               %    mfit=array(number of objectives by 
                               %               number. of ind. in region)
                               %               of raw fitness values
                               %    fit=vector with aggregate fitness 
                               %        values of individuals in region                  
                  
% Death algorithm parameters
GAP.dt_alg=2;                  % algorithm
                               %    1 = replace parents
                               %    2 = random replacement
                               %    3 = tournament on fitness
                               %    4 = tournament on age
                               %    5 = custom
                               %    6 = random (select 1-4)
GAP.dt_nts=4;                  % number of individuals in death 
                               %    tournament selection    
GAP.dt_cah=[];                 % custum algorithm handle; should be of form
                               % death_list=f(region,size,age,mfit,fit)
                               % where
                               %    death_list=indices of those individuals 
                               %               to be replaced by children
                               %    region=region number
                               %    size=number of individuals to be in 
                               %         death list
                               %    age=vector describing ages of 
                               %        individuals in population
                               %    mfit=array (number of objectives by 
                               %         number of individuals in region)
                               %         with raw fitness values
                               %    fit=vector with aggregate fitness
                               %         values of individuals in region

% Mating and crossover parameters
GAP.mc_pp   = 0.6;             % percentage of population replaced 
                               %    by children
GAP.mc_fc   = 1.0;             % probability (fraction) of chromosomes of 
                               %    an individual who is selected 
                               %    for mating to undergo crossover
GAP.mc_alg  = 4.0;             % crossover algorithm
                               %    1 = single point crossover
                               %    2 = scalar simple blend crossover
                               %    3 = vector simple blend crossover
                               %    4 = scalar simulated binary crossover
                               %    5 = vector simulated binary crossover
                               %    6 = random algorithm
if (GAP.fp_obj~=0)
   GAP.mc_alg=6;
end
GAP.mc_gac  = 3.0;             % generations between algorithm 
                               % change for random crossover
GAP.mc_ec   = 2.0;             % Distribution tightness parameter eta sub c 
                               % for simulated binary crossover

% Parameters associated with mutation
GAP.mt_ptgm0=0.01;   % initial probability of total gene mutated
GAP.mt_ptgmf=0.001;  % final probability of total gene mutated
GAP.mt_prgm0=0.02;   % init. probability of relative partial gene mutation
GAP.mt_prgmf=0.002;  % fin. probability of relative partial gene mutation
GAP.mt_srgm0=0.3;    % int. std. deviation of relative partial gene mutation
GAP.mt_srgmf=0.03;   % fin. std. deviation of relative partial gene mutation
GAP.mt_pagm0=0.02;   % init. prob. of absolute partial gene mutation
GAP.mt_pagmf=0.002;  % fin. prob. of absolute partial gene mutation
GAP.mt_sagm0=0.1;    % init. std. dev. of absolute partial gene mutation
GAP.mt_sagmf=0.01;   % fin. std. dev. of absolute partial gene mutation
GAP.mt_prvm0=0.02;   % init. prob. of relative vector mut. of individual
GAP.mt_prvmf=0.002;  % fin. prob. of relative vector mutation of individual
GAP.mt_srvm0=0.3;    % init. std. dev. of relative vector mutation
GAP.mt_srvmf=0.03;   % fin. std. dev. of relative vector mutation
GAP.mt_pavm0=0.02;   % init. prob. of absolute vector mut. of individual
GAP.mt_pavmf=0.002;  % fin. prob. of absolute vector mutation of individual
GAP.mt_savm0=0.1;    % init. std. dev. of absolute vector mutation
GAP.mt_savmf=0.01;   % fin. std. dev. of absolute vector mutation
GAP.mt_pigm0=0.2;    % init. prob. of integer gene mutation
GAP.mt_pigmf=0.01;   % fin. prob. of integer gene mutation

% Gene repair
GAP.gr_alg=1;     % algorithm
                  % 1 = hard limit
                  % 2 = ring mapping

% Migration parameters
GAP.mg_nreg=round(GAP.fp_npop/100);  % number of regions
if GAP.mg_nreg<1
   GAP.mg_nreg=1;
end
GAP.mg_tmig=round(GAP.fp_ngen/20);   % time between migrations in gen.
if (GAP.mg_tmig<1)
   GAP.mg_tmig=1;
end
GAP.mg_pmig=0.1;                     % probability of individual 
                                     % migrating when opportunity arises

% Evaluation parameters                               
GAP.ev_bev=0;     % if 0 single evaluation; if 1 block evaluation
GAP.ev_are=0;     % if 0: evaluate only the unevaluated individuals
                  % if 1: evaluate all individuals
GAP.ev_ssd=0;     % send supplementary data
                  %    if 0: don't send anything to fitness function 
                  %          other than genevalues and, optionally,
                  %          externally supplied data
                  %    if 1: send region number and previous fitness 
                  %          raw fitness value
GAP.ev_pp=false;  % parallel process
GAP.ev_npg=12;    % number of evaluation groups for non-block
                  % evaluation parallel processing
                                               
% Elitism parameters
GAP.el_act=1;     % activate elitism 
GAP.el_fgs=0.0;   % fraction of generations through study to start elitism
GAP.el_fpe0=0.2;  % fraction of population protected as elite for 
GAP.el_fpef=0.8;  % multi-objective optimization
                        
                           
% Random search parameters
GAP.rs_fgs=0.5;   % fraction of gen. through study to start random search
GAP.rs_fps=0.1;   % fraction of total population size used in random search 
GAP.rs_srp=0.3;   % standard deviation used in relative perturbation
GAP.rs_sap=0.1;   % standard deviation used in absolute perturbation
GAP.rs_frp=0.7;   % fraction of relative random perturbations
                  % fract. of abs. random perturbations is (1-GAP.rs_frp)
GAP.rs_fea=0.2;   % fract. of geneations on which to execute the algorithm

% Reporting parameters
GAP.rp_lvl=1;                  % reporting level 
                               % -1 = no plots or reports
                               %  0 = text report only
                               %  1 = plots and text reporting
GAP.rp_gbr=5;                  % generations between reports
GAP.rp_crh=[];                 % handle to custom reporting function
                               % example use:  GAP.rp_crh=@myreport
                               % where myreport(P,GAP) is where the user 
                               % function would be called
                               % and where P is the current population.  
                               % Note this will be
                               % called every GAP.rp_gbr generations                               
                               
% Objective plot parameters
if GAP.fp_obj~=0
    GAP.op_list=[GAP.fp_obj];  % list of objectives to make obj. plots for
else
    GAP.op_list=[];
end 
GAP.op_style(1:GAP.fp_nobj)=1; % style for each objective (0=log, 1=linear)
GAP.op_sign(1:GAP.fp_nobj)=1;  % sign of fitness for each objective 
                               % (1=positive/mixed,-1=negative)

% Parameter to determine which obj. the output plot should be sorted by.
GAP.op_srt = 1;                % Default is the first objective. 


% Distribution plot parameters
GAP.dp_type=1;                 % type of distribution plot
                               %    1 = plot individuals
                               %    2 = plot histograms
GAP.dp_np=100;                 % max. number of pop. to plot for type 1
GAP.dp_res=20;                 % number of bins in distr. plot for type 2

% Pareto plot parameters
if GAP.fp_obj==0
   GAP.pp_list=[1,2];          % list of 2 or 3 pars. used in Pareto plot
else
   GAP.pp_list=[];
end
GAP.pp_xl='Objective 1';       % x-axis label
GAP.pp_yl='Objective 2';       % y-axis label
GAP.pp_zl='Objective 3';       % z-axis label
GAP.pp_title='Solution Space'; % Pareto plot title
GAP.pp_style(1:GAP.fp_nobj)=1; % style for each obj. (0=log, 1=linear)
GAP.pp_sign(1:GAP.fp_nobj)=1;  % sign of fitness for each objective 
                               % (1=positive/mixed,-1=negative)
GAP.pp_axis=[];                % axis limits for Pareto plot

% Output parameters 
GAP.so_ob = 1; % This points to the objective to be used to sort the output
               % GAP.so_ob = 1  - Sorts with respect to Objective 1. 
       
% trimga
GAP.trimga=1;                  % trim ga results on single objective

% The following fields must be defined by user
% GAP.gd_min                   % row vector of minimum gene values
% GAP.gd_max                   % row vector of maximum gene values
% GAP.gd_type                  % row vector of gene types
                               % 1 = integer
                               % 2 = real
                               % 3 = logarithmic
% GAP.gd_cid                   % row vector of chromosome id

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