function [fP,GAS,bI,f]= gaoptimize(objhandle,GAPic,D,GAS,iP)
% GAOPTIMIZE  Optimizes a function using a genetic algorithm.
%
% [fP,GAS]      = gaoptimize(@fitfun,GAP,D,GAS,iP)
% [fP,GAS]      = gaoptimize(@fitfun,GAP,D)
% [fP,GAS]      = gaoptimize(@fitfun,GAP)
% [fP,GAS,bI]   = gaoptimize(@fitnun,GAP,D,GAS,iP)
% [fP,GAS,bI]   = gaoptimize(@fitnun,GAP,D)
% [fP,GAS,bI]   = gaoptimize(@fitnun,GAP)
% [fP,GAS,bI,f] = gaoptimize(@fitnun,GAP,D,GAS,iP)
% [fP,GAS,bI,f] = gaoptimize(@fitnun,GAP,D)
% [fP,GAS,bI,f] = gaoptimize(@fitnun,GAP)
%
% Inputs:
% fitfun   = name of *.m function that evaluates the fitness.  If blckeval 
%            is set to zero (see below) the input is a vector of length 
%            equal to the number of genes, and the output is a column 
%            vector of fitness functions of length equal to the number of 
%            objectives. If blckeval is set to 1 for block evaluation the 
%            input is a matrix where the each row corresponds to a 
%            different gene and each column is a different individual.  In 
%            this case the output is also a matrix where each row 
%            corresponds to a different objective and each column is a 
%            different member of the population.
% GAP      = Genetic Algorithm Parameters
%            (see gap_default for definitions)
% D        = data needed by fitness function (could be an empty matrix)
% GAS      = Genetic Algorithm Statistics (use empty matrix when not used)
% iP       = Optional variable with intial population (Use empty matrix 
%            when not in use)
%
% Outputs:
% bI       = Best individual (For single objective case, these are the gene
%            values for the best individual. For the multi objective case,
%            each column is a non-dominated individual (sorted according to
%            the specified objective).
% fP       = Final population
% GAS      = Genetic Algorithm Statistics
% f        = Best fitness/fitnesses(Each column is a set of fitness values
%            of the individuals corresponding to bI)
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

% alternate input form for GAP
if isfield(GAPic,'gd')
    GAPic.gd_min=GAPic.gd(:,1)';
    GAPic.gd_max=GAPic.gd(:,2)';
    GAPic.gd_type=GAPic.gd(:,3)';
    GAPic.gd_cid=GAPic.gd(:,4)';
end

% initialize GAP
GAP=GAPic;

% set up
switch nargin
    case 5
    case 3
        GAS=[];
        iP=[];
    case 2
        D=[];
        GAS=[];
        iP=[];
    otherwise
        error('Incorrect number of arguments to gaoptimize');
end

switch nargout
    case 2
        bI=[];
    case 3
        f=[];
    case 4
    otherwise
        error('Incorrect number of output arguments to gaoptimize');
end

% inform the user as to the study we are going into
if (GAP.fp_obj~=0)
   disp(['Performing single-objective optimization of objective ' ...
       num2str(GAP.fp_obj) ' using a population size of ' ...
       num2str(GAP.fp_npop) ' over ' num2str(GAP.fp_ngen) ' generations']);
else
   disp(['Performing multi-objective optimization ' ...
         'using a population size of ' num2str(GAP.fp_npop) ' over ' ...
         num2str(GAP.fp_ngen) ' generations']);
end

% initialization
[GAP,GAS,P] = gainit(objhandle,D,GAP,GAS,iP);
initgen = GAS.cg;

% let the population evolve
for j=initgen:GAP.fp_ngen,
    
    % Update current population number
    GAS.cg = j;                       
    
    % Update GA parameters
    GAP=gapadjust(GAPic,GAS);        
    
    % Objective weighting vector for current generation
    GAP.owv                     = objwght(GAP);    
    
    % Diversity dontrol
    P.pen                       = divcon(P,GAP); 
    
    % Determine the scaled fitness
    P.fit                       = scale(P,GAP); 
    
    % Select mating pool / parent list
    [P_list,PL_size]            = select(P,GAP); 
    
    % Select individuals to die
    D_list                      = death(P,P_list,GAP);   
    
    % Mating and Crossover
    P2                          = matingcrossover(P, ...          
                                  P_list,PL_size,D_list,GAP,GAS);
                              
    % Mutation                          
    P3                          = mutate(P2,GAP);         
    
    % Migration
    P4                          = migrate(P3,GAP,GAS.cg);        
    
    % Update age
    P4.age                      = updateage(P4);        
    
    % Evaluate
    [P4.mfit,P4.eval,GAS.ne(j)] = evaluate(P4,GAP,GAS.ne(j-1),D); 
    
    % Elitism
    P5                          = elitism(P4,P,GAP,GAS);    
    
    % Random Search
    [P,GAS]                     = randsearch(P5,GAP,GAS,D);    
    
    % Update statistics
    [GAS] = updatestat(GAP,GAS,P); 
    
    % Update report plot
    reportplot(GAP,GAS,P);                                 
   
end

fP = P;      % Final population



% Updated 1 July 07
% Describe simulation
if (GAP.fp_obj~=0)
    disp(['This single-objective optimization of objective ' ...
        num2str(GAP.fp_obj) ' using a population size of ' ...
        num2str(GAP.fp_npop) ' over ' num2str(GAP.fp_ngen) ...
        ' generations']);
else
    disp(['This multi-objective optimization ' ...
        'used a population size of ' num2str(GAP.fp_npop) ' over ' ...
        num2str(GAP.fp_ngen) ' generations']);
end

disp(' ');

% Updated 1 July 07
% Pick out and return best individual
if (nargout==3)||(nargout==4)
    
    if (GAP.fp_obj>0)
    
        % Single-Objective Case
        if GAP.trimga
            disp('Trimming best answer ...');
            if ~isempty(D)            
                [bI,f]=trimga(GAP,fP,D);
            else
                [bI,f]=trimga(GAP,fP);
            end
            bI=bI';
        else
            bI=GAS.bestgenes(:,end,GAP.fp_obj);
        end
        
    else
        
        % Multi-Objective Case
        % Identify non-dominated solutions
        ind = logical(nondom(P.mfit,1)); 
        fit_nd = P.mfit(:,ind);
        gene_val = P.gene(:,ind);
       
        % GAP.fp_nobj is the number of objectives.
        if GAP.fp_nobj < GAP.so_ob
            error('Wrong input to GAP.so_ob')
        end

        % Sorting based on user input.
        [~,ind] = sort(fit_nd(GAP.so_ob,:));
        bI = gene_val(:,ind); % bI is a matrix with non-dominated 
                              %solutions' gene values
        f = fit_nd(:,ind);    % f is the objective values
        
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