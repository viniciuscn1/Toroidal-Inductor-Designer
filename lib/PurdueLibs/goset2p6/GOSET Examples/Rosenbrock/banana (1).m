% Optimization of Rosenbrock's Problem

% Initialize the parameters
GAP = gapdefault;                            % Line 1

% Define gene parameters 
%                 x1     x2    
% gene             1      2     
GAP.gd_min  = [   -2     -1    ];            % Line 2
GAP.gd_max  = [    2      3    ];            % Line 3
GAP.gd_type = [    2      2    ];            % Line 4
GAP.gd_cid  = [    1      1    ];            % Line 5

% Execute GOSET
[P,GAS,best] = gaoptimize(@bananafit,GAP);   % Line 6
