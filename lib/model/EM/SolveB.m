function B = SolveB(P,Ic,r)
%%
% This function compute flux density across the toroidal inductor with
% nonlinear anhysteretic BH curve, and permeability multiplier.
% 
% Call: B = SolveB(P,Ip,r)
% 
% Input:
%  P        = parameters structure
%   P.mur   = permeability multiplier piece function
%   P.N     = # of turns in the toroid coil
%  Ic       = coil current [A]
%  r        = radii test points [m]
% 
% Output:
%  B        = flux density for every radii point [T]
% 
% Internal variables:
%  k        = permeability multiplier
%  B0       = initial guess of flux density for the nonlinear solve [T]
% 
% Written by Vinicius Cabral do Nascimento
%   email: vcabrald@purdue.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve nonlinear core domain
r = r(:);                       % enforce column vector for radius [m]
Ic= transpose(Ic(:));           % set current as a row vector [A]
k = ppval(P.mur,r);             % compute permeability multiplier
B0= zeros(length(r),length(Ic));% initialize B0
B = fsolve(@(B) ConstitutiveResidual(P,Ic,r,k,B),B0,optimoptions(...
    'fsolve','SpecifyObjectiveGradient',true,'Display','none'));

end

%%
% private functions
function [F,dF] = ConstitutiveResidual(P,Ip,r,k,B)
%%
% This function computes the magnetic constitutive relations residual
% 
% Call: [F,dF] = ConstitutiveResidual(P,Ip,r,k,B)
% 
% Input:
% 
% Output:
% 
% Internal variables:
%  mu       = magnetic permeability [H/m]
%  pmu      = derivative of magnetic permeability w.r.t. flux density
%  H        = field intensity [A/m]
%  F        = residual function of B-mu*H
%  dF       = derivative of residual function w.r.t. flux density
%
% Written by Vinicius Cabral do Nascimento
%   email: vcabrald@purdue.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% get mu(B) and pmu(B)
[mu,pmu] = fmuB(P.mp,B,k);

% compute field intensity
H   = (P.N*Ip(ones(length(r),1),:))./(2*pi*r(:,ones(1,length(Ip))));
F   = B-mu.*H;          % get function residual

% in case Jacobian matrix is requested
if nargout>1
    dF = 1-pmu.*H;
    dF = diag(dF(:));
end

end % end of ConstitutiveResidual(P,Ip,r,k,B)