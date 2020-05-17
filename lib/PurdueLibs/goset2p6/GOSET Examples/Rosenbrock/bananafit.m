function f = bananafit(x)
% Rosenbrock's banana function 

% Banana function
f1 = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;  

% Fitness function
f = 1/(0.001 + f1);                         