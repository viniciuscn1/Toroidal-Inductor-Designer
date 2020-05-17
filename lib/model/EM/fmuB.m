%  FILE:          muB.m
%  FUNCTION:      [mu,pmu] = muB(MP,B)
%  DESCRIPTION:   This routine calculates permeability as a function of B.
%
%  Inputs:       MP       = Structure of material parameters
%                 MP.mur  = Initial relative permeability
%                 MP.muB.a[]  = Vector of alpha coefficients (1/T)
%                 MP.muB.b[]  = Vector of beta exponential coefficients
%                               (1/T)
%                 MP.muB.g[]  = Vector of gamma exponential offsets (T)
%                 MP.muB.d[]  = Vector of delta coefficients
%                 MP.muB.e[]  = Vector of epsilon values
%                 MP.muB.z[]  = Vector of zeta values
%                 MP.muB.h[]  = Vector of eta values (1/T)
%                 MP.muB.t[]  = Vector of theta values
%                B        = Vector of points at which to calculate 
%                           permeability (T)
%
%  Outputs:      mu       = Permeability at points corresponding to B (H/m)
%                pmu      = derivative of mu with respect to B (H/(m*T))
%
% Internal:      aB       =  absolute value of B (T)
%                mu0      =  permeability of freespace (H/m)
%                NT       =  number of terms
%                f        =  the ratio of flux density over magnetization 
%                            (with both expressd in T)
%                pg       =  derivative of g with respect to B (1/T)
%                            (g=f less a constant)
%                n        =  index variable
%                Bterm    =  term in expression
%
%  AUTHOR:      Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%  REFERENCE: G.M. Shane, S.D. Sudhoff, “Refinements in Anhysteretic
%             Characterization and Permeability Modeling,”  IEEE 
%             Transactions on Magnetics, vol. 46, no. 11, pg(s) 3834-3843, 
%             November 2010.
function [mu,pmu] = fmuB(MP,B,k)

%%
% process input
k = k(:);               % ensure permeability tuning as column vector

% define local variales
aB  = abs(B);           % absolute value of the flux density
mu0 = 4e-7*pi;          % magnetic permeability of free space [H/m]
nt  = length(MP.muB.a); % # of terms in the permeability function

% initialize f and pg
f   = k*MP.mur./(k*MP.mur-1);   
pg  = 0;

% build f and pg
for n=1:nt              % loop through each term of the permeability function
   Bterm= exp(-MP.muB.b(n)*aB); 
   f    = f+MP.muB.a(n)*aB+MP.muB.d(n)*log(MP.muB.e(n)+MP.muB.z(n)*Bterm);
   pg   = pg+MP.muB.h(n)./(MP.muB.t(n)+Bterm);
end

% define output
mu  = mu0*f./(f-1);     % permeability at given magnetization
pmu = -mu0*sign(B).*pg./(f-1).^2;   % derivative of magnetic permeability
end