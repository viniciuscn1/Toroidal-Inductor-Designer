% Transfer Function Identification Example
%
% 	Stand Still Frequency Response (SSFR) 
%   Impedance Data of Brushless DC Machine
% 	Rotor is aligned with d-axis.  

% Get admittance data
ydata;        
data.s = s;
data.t = yd; 

% Initialize the gemetic algorithm parameters
GAP = gapdefault(1,1,200,400);  % Set default values for GAP
GAP.mc_alg = 6;
GAP.dt_alg = 3;

% Set range for genes
order = 6;
O=ones(1,order);
GAP.gd_min =  [ 1e-8*O 1e-8*O ];
GAP.gd_max =  [ 1e+1*O 1e+0*O ];
GAP.gd_type = [ 3*O 3*O ];
GAP.gd_cid  = [ 1*O 1*O ];

% Execute GOSET
[P,GAS, best] = gaoptimize(@tfid_fit,GAP,data);

% Plot the results
a   = best(1:order);
tau = best(order+1:2*order);
ydp = pftf(s,a,tau);
figure(2);
bodeplot(2,f,yd,'bx',ydp,'r-');


