function pld = ComputeMSE(qe,B,we,CM)
% ComputeMSE computes the core loss density of a material
%            based on a 1/4 cycle data. Assumes B is symmetric
%            and half-wave symmetric. Uses a variation of the
%            modified Steinmetz equation.
%
% Call:
% pld=core_loss_density(qe,B,we,CM)
%
% Inputs:
% qe         = angle of data ranging over 1/2 of a cycle (rad)
% B          = flux density structure with data over 1/2 of a cycle (T)
% we         = radian frequency of data
% CM         = stucture of material data
% 
% Outputs:
% Pld        = power loss density (W/m^3)
% 
% Internal:
% f          = fundamental frequency (Hz)
% t          = time (s) (spans 1/4 cycle)
% N          = number of points in vectors
% i1,i2      = index arrays
% pB         = time derivative of B (T/s)
% pB2        = square of time derivative of B (T/s)^2
% Bmx        = maximum flux density (T)
% Bpk        = peak flux density (T)
% DB         = peak-to-peak flux density range (T)
% int_dBdt2  = integral of square of derivative of B over a cycle
% feq        = equivalent frequency
%
% Written by:
% Vinicius Nascimento
% S.D. Sudhoff

% define local constants
np = length(qe);            % # of points in the waveform
f = abs(we)/(2*pi);         % fundamental frequency [Hz]
t = abs(qe/we);             % elapsed time [s]
t = t-t(1);                 % adjust beginning of time vector to zero 

% compute time derivative of flux density
pB = diff(B,1,2)./repmat(diff(t),[size(B,1),1]);
pB(:,np) = pB(:,np-1);
   
% Determine the square of the time derivative of flux density
pB2 = pB.^2;                     
   
% Determine maximum B
Bmx = max(abs(B),[],2);
Bpk = Bmx;
DB  = 2*Bmx;
    
% Determine integral (factor of 4 is because using 1/4 of a cycle)
int_dBdt2 = 4*trapz(t,pB2,2);
  
% Determine equivalent frequency
feq = zeros(size(DB));
feq(DB~=0) = 2*int_dBdt2(DB~=0)./(DB(DB~=0)*pi).^2;
               
% Determine output
pld = CM.SE.kh*Bpk.^(CM.SE.beta).*feq.^(CM.SE.alpha-1)*f + ...
    CM.SE.ke*f*int_dBdt2;

end