% Tanaka problem (1995)
%
% This is a multi-objective test function, with unconstrained domain;
% proposed by Poloni (MOP3, p.110 [1])
%
% Reference:
% [1] C. A. Coello Coello, D. A. Van Veldhuizen, and G. B. Lamont,
%     "Evolutionary Algorithms for Solving Multi-Objective Problems"
%     Kluwer Academic Publishers, 2002
%
% date: 12.5.2003

function [f] = tanaka_fit(x)

C1 = x(1)^2+x(2)^2-1-0.1*cos(16*atan(x(1)/x(2))) >= 0;
C2 = (x(1)-0.5)^2+(x(2)-0.5)^2 <= 0.5;

if C1 && C2   
    f(1,1) = -x(1);
    f(2,1) = -x(2);
else
    f(1,1) = -10;
    f(2,1) = -10;
end
