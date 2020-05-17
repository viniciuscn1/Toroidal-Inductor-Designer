% This is a multi-objective test function, with unconstrained domain;
% proposed by Poloni (MOP3, p.110 [1])
%
% Reference:
% [1] C. A. Coello Coello, D. A. Van Veldhuizen, and G. B. Lamont,
%     "Evolutionary Algorithms for Solving Multi-Objective Problems"
%     Kluwer Academic Publishers, 2002
%
% date: 12.5.2003

% Initialize the parameters
GAP = gapdefault(2,0,200,200);

% Plotting parameters
GAP.op_list = [];               % objectives list for objective plots
GAP.pp_list = [ 1, 2];          % gene list for Pareto plot
GAP.pp_sign = [-1,-1];          % sign of fitness for each objective 
                                % (1=positive,-1=negative)
GAP.pp_axis = [0 1.25 0 1.25];  % axis limits for Pareto plot

% Gene setup
%                      x1          x2    
% gene                  1           2     
GAP.gd_min  = [    0.0001      0.0001     ];
GAP.gd_max  = [        pi          pi     ];
GAP.gd_type = [         2           2     ];
GAP.gd_cid  = [         1           1     ];
  
% Perform optimization
[P,GAS,best]= gaoptimize(@tanaka_fit,GAP);

% Plot final solutions
figure(2)
plot(best(1,:),best(2,:),'x');
axis([0 1.25 0 1.25]);
xlabel('f_1');
ylabel('f_2');
axis square;
title('Pareto Front')