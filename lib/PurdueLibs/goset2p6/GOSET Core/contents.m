% Genetic Optimization and System Engineering Tool (GOSET)
% Version: 2.6
%
% S.D. Sudhoff
% Purdue University
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% E-mail: sudhoff@ecn.purdue.edu
%
% Common Data Stuctures 
%
% GOSET 2.6 is built around 3 data structures - P, GAP, and GAS.
% P represents a population of chromosomes.  GAP is a structure
% containing information on the parameters of the genetic algorithm.
% GAS is a structure containg information on the statistics of the
% genetic algorithm.
%
% P is the structure for a population. Its fields consist of 
%
% P.blckeval   If this exist, and set to 1, then the fitness function 
%              expects its argument to include a set of studies rather
%              than a single study.
% P.fithandle  Handle to function to be used to evaluate the population
% P.size       Size (number of chromosomes) in population.
% P.nobj       The number of objective functions to maximize.
% P.mfit       An array P.nobjxP.size desscribing the multi-objective 
%              fitness of each member in the population
% P.fit        An vector 1xP.size which is a scaled, penalized, and 
% P.eval       An array 1xP.size describing weather the fitness
%              aggegated fitness of that member is current (1) or if it 
%              needs to be re-evaluated (0).
% P.age        A vector 1xP.size which describes the age of the population
% P.ngenes     Number of genes.
% P.min        Array (P.ngenes x 1) of minimum value of gene
% P.max        Array (P.ngenes x 1) of maximum value of gene
% P.type       Array (P.ngenes x 1) of type of genes 
%              (1=integer,2=real,3=exponential)
% P.chrom_id   Array (P.ngenes x 1) describing which chromosome a gene is on
% P.normgene   Array (P.ngenes x P.size) of normalized gene values
% P.gene       Array (P.ngenes x P.size) of gene values
% P.region     Geographic region (1 x P.size) of a population member
% P.pen        An array 1xP.size describing the penalty function of each member in
%              the population
%
% GAP is a structure for the parameters of the genetic algorithm.  It is
% documented in the gap_default routine which sets default values.
%
% GAS is a structure with the genetic algorithm statistics.  Its fields are
% as follows:
%
% GAS.cg          Current generation
% GAS.medianfit   A matrix (number of objectives x number of generations) 
%                 describing the median fitness of each objective
% GAS.meanfit     A matrix (number of objectives x number of generations) 
%                 describing the mean (average) fitness of each objective
% GAS.bestfit     A matrix (number of objectives x number of generations) 
%                 describing the best fitness of each objective
% GAS.bestgenes   A three-diminsional matrix (number of genes x number of 
%                 generations x number of objectives) whose descring the 
%                 best gene values for generation for each objective
% GAS.ne          The number of objective function evaluations made
%
% GOSET consists of a number of routines to support the use of genetic
% algorithms to solve problems.  The most commonly used routines are:
%
% gapdefault   This routine sets up the default values for all the
%              parameters of the genetic algorithm
% gaoptimize   This routine performs a genetic optimization of an arbitrary
%              function
%
% A complete list of the funtions in GOSET
% include
%
% distplot        Plots gene distribution.
% divcon          Diversity control.
% death           Determines individuals to die (replaced by parents)
% downsize        Reduces the size of a population.
% elitism         Used to implement the elitism operater.
% goset           Graphical user interface for the toolbox.
% evaluate        Evaluates the fitness funtion of a pop. of chromsomes.
% gapadjust       Determines mutation parameters for current generation
% gainit          Performs initialization functions. 
% gaoptimize      Optimizes and objective function using a GA.
% gapdefault      Sets up default parameters of the genetic algorithm used.
% gasetup         Sets up a data structure P for a pop. of chromosomes.
% gosetinit       Initializes graphical user interface.
% matingcrossover Mating and crossover.
% migrate         Performs migration operation between regions.
% mutate          Mutation operator.
% nondom          Find non-dominated solutions.
% normgene        Updates normalized gene values based on raw gene values.
% objwght         Objective weighting function.
% paretoplot      Plots the individuals in objective function space so the
%                 evolution of the pareto optimal front can be observed.
% randsearch      Performs a random search about most fit individual in 
%                 each region.
% rawgene         Updates raw gene values based on normalized genes.
% reportplot      A routine for reporting and plotting progress of a
%                 genetic algorithm.
% scale           Scales, penalizes, and aggregates fitness.
% select          Selects a parent list
% unrndinit       Uniform random initialization of a population.
% updateage       Updates age of a population
% updatestat      Update statistics.

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