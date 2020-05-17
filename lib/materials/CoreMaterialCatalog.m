function [SP] = CoreMaterialCatalog(n,fcl)

% CoreMaterialCatalog assigns the core parameters to a parameter structure
%
% [SP] = CoreMaterialCatalog(n,fcl)
%
% Inputs
% n             = material index 
%                 n=1 is Strain Anneal Fe=Mn=2.7
%                 n=2 is Field Anneal Co-Based 300perm
%                 n=3 is Powder Core MPP14
% fcl           = frequency to compute core losses [Hz]
% 
% Outputs
% SP            = structure of core material parameters 
%  SP.desc      = core material description
%  SP.Msat      = saturated magnetization (T)
%  SP.mur       = relative permeability in linear region
%                 (from anhysteretic results)
%  SP.Blim      = recommended limit on B to avoid saturation (T) 
%                 This is taken the point where the absolute relative 
%                 permeability hits mu_int/100.
%  SP.rho       = mass density [kg/m^3]
%  SP.BH        = structure of BH curve parameters
%   SP.BH.m     = magnetization coefficients (T)
%   SP.BH.n     = exponents
%   SP.BH.h     = field intensity breakpoints (A/m)
%  SP.MSE       = structure of MSE loss parameters
%   SP.SE.alpha= frequency exponent 
%   SP.SE.beta = flux density exponent
%   SP.SE.kh   = hysteresis loss coefficient (J/m^3)
%   SP.SE.ke   = eddy current loss coefficient (J*s/m^3)
%  SP.muB       = structure of mu(B) parameters
%   SP.muB.mur  = initial relative permeability of anhysteric curve                                anhysteretic curve
%   SP.muB.a    = vector of alpha coefficients (1/T)
%   SP.muB.b    = vector of beta exponential coefficients (1/T)
%   SP.muB.g    = vector of gamma exponential offsets (T)
%   SP.muB.d    = vector of delta coefficients
%   SP.muB.e    = vector of epsilon values
%   SP.muB.z    = vector of zeta values
%   SP.muB.h    = vector of eta values (1/T)
%   SP.muB.t    = vector of theta values
%  SP.TEC       = thermal properties structure
%   SP.TEC.kx   = thermal conductivity in x [W/(m*K)]
%   SP.TEC.ky   = thermal conductivity in y [W/(m*K)]
%   SP.TEC.kz   = thermal conductivity in z [W/(m*K)]
%   SP.TEC.c    = specific heat capacity [J/(kg*K)]
%
% References:
% Shane, G. and S. D. Sudhoff,  "Refinements in Anhysteretic
% Characterization and Permeability Modeling,"  IEEE
% Transactions on Magnetics, vol. 46, no. 11 November 2010.

SP.pt_muB = @muB;

if nargin<2, fcl = 0; end

switch n
    
    case 1 % 'Strain Anneal Fe=Mn=2.7' ------------------------------------
        
        % general parameters
        SP.desc = 'Strain Anneal Fe=Mn=2.7';
        SP.Msat = 1.3689;       % saturation magnetization
        SP.mur  = 30.5658;      % relative permeability
        SP.Blim = 1.933;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters as a function of H
        SP.BH.m = [1.773,       0.72031,    0.2397,    0.31238];
        SP.BH.n = [1,           3.1215,     3.3259,    5.9978];
        SP.BH.h = [118129.5423, 53013.5513, 125546.72, 34449.62345];
        
        % MSE parameters
        SP.SE.alpha= 1.08102308240564;
        SP.SE.beta = 2.01596367003735;
        SP.SE.kh   = 45.616848607384;
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 30.5658;
        SP.muB.a = [0.3396    0.023712   0.0063699   0.0015238];
        SP.muB.b = [56.03812      61.81279      26.04038          1000];
        SP.muB.g = [1.0284     0.86461      1.1314     0.46935];
    
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 7;      % thermal conductivity in z [W/(m.K)]
        
    case 2 % 'Field Anneal Co-Based 300perm' ------------------------------
        
        % general parameters
        SP.desc = 'Field Anneal Co-Based 300perm';
%         SP.Msat = 1.8571;
        SP.mur  = 350;
%         SP.Blim = 1.4874;
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.8571     0.83433    -0.51725      1.1778];
%         SP.BH.n = [1       3.304      1.2606      2.0209];
%         SP.BH.h = [289.999       159.604      160.2678      296.7384];
        
        % MSE parameters
        SP.SE.alpha= 1.053995;
        SP.SE.beta = 2.094071;
        SP.SE.kh   = 10.06432; 
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 350;
        SP.muB.a   = 0.55;
        SP.muB.b   = 35.0;
        SP.muB.g   = 0.96;

        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 7;      % thermal conductivity in z [W/(m.K)]

    case 3 % 'Powder Core MPP14' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP14';
        SP.Msat = 1.4471;
        SP.mur  = 14;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 64.02;
            b = 1.074;
            c = 1.11;
        else
            a = 21.06;
            b = 1.074;
            c = 1.38;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1); 
        SP.SE.ke   = 0;

        % mu(B) parameters
        SP.muB.mur = 14;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
        
    case 4 % 'Powder Core MPP26' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP26';
        SP.Msat = 1.4471;
        SP.mur  = 25.9335;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 361.62;
            b = 2.000;
            c = 1.08;
        else
            a = 109.17;
            b = 2.000;
            c = 1.37;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1); 
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 25.8616;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
        
    case 5 % 'Powder Core MPP60' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP60';
        SP.Msat = 1.4471;
        SP.mur  = 60;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 80.12;
            b = 1.585;
            c = 1.04;
        else
            a = 31.32;
            b = 1.585;
            c = 1.37;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1); 
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 60;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 6 % 'Powder Core MPP125' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP125';
        SP.Msat = 1.4471;
        SP.mur  = 125;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 254.26;
            b = 2.222;
            c = 1.17;
        else
            a = 87.07;
            b = 2.222;
            c = 1.56;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 125;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 7 % 'Powder Core MPP147' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP147';
        SP.Msat = 1.4471;
        SP.mur  = 147;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 254.26;
            b = 2.222;
            c = 1.17;
        else
            a = 87.07;
            b = 2.222;
            c = 1.56;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 147;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 8 % 'Powder Core MPP160' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP160';
        SP.Msat = 1.4471;
        SP.mur  = 160;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 254.26;
            b = 2.222;
            c = 1.17;
        else
            a = 87.07;
            b = 2.222;
            c = 1.56;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 160;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 9 % 'Powder Core MPP173' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP173';
        SP.Msat = 1.4471;
        SP.mur  = 173;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 254.26;
            b = 2.222;
            c = 1.17;
        else
            a = 87.07;
            b = 2.222;
            c = 1.56;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 173;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 10 % 'Powder Core MPP200' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP200';
        SP.Msat = 1.4471;
        SP.mur  = 200;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 320.32;
            b = 2.322;
            c = 1.19;
        else
            a = 115.52;
            b = 2.322;
            c = 1.59;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 200;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 11 % 'Powder Core MPP300' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP300';
        SP.Msat = 1.4471;
        SP.mur  = 300;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 320.32;
            b = 2.322;
            c = 1.19;
        else
            a = 115.52;
            b = 2.322;
            c = 1.59;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
        
        % mu(B) parameters
        SP.muB.mur = 300;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    case 12 % 'Powder Core MPP500' ------------------------------------------
        
        % general parameters
        SP.desc = 'Powder Core MPP500';
        SP.Msat = 1.4471;
        SP.mur  = 500;
        SP.Blim = 1.887;        % B corresponding to mur_init/10
        SP.rho  = 5500;         % density [kg/m3]
        
        % BH parameters
%         SP.BH.m = [1.4471  -0.0019674     0.47977    -0.65306];
%         SP.BH.n = [1      1.8908      2.2776      1.2747];
%         SP.BH.h = [24769.1645      1937.04522      26969.9187      15326.6643];
        
        % MSE parameters: "1" is when frequency <= 10kHz, "2" is when freq > 10kHz
        if fcl<=10e3
            a = 303.43;
            b = 1.999;
            c = 1.09;
        else
            a = 96.89;
            b = 1.999;
            c = 1.54;
        end
        SP.SE.alpha= c;
        SP.SE.beta = b;
        SP.SE.kh   = a/1e3^(c-1);
        SP.SE.ke   = 0;
                
        % mu(B) parameters
        SP.muB.mur = 500;
        SP.muB.a = [0.085658    0.012158    0.011584    0.005437];
        SP.muB.b = [21.60639          1000      70.66507          1000];
        SP.muB.g = [0.5983     0.26256     0.43251     0.34387];
        
        % TEC parameters
        SP.TEC.kx = 1.5;    % thermal conductivity in x [W/(m.K)]
        SP.TEC.ky = 1.5;    % thermal conductivity in y [W/(m.K)]
        SP.TEC.kz = 1.5;    % thermal conductivity in z [W/(m.K)]
    otherwise %------------------------------------------------------------
        
        % assign NaN to all parameters and exit function
        SP.desc='BAD';
        SP.Msat=NaN;
        SP.mur=NaN;
        SP.Blim=NaN;
        SP.row=NaN;
        SP.k=NaN;
        SP.c=NaN;
        
        SP.BH.m=NaN;
        SP.BH.h=NaN;
        SP.BH.n=NaN;
        
        SP.SE.alpha=NaN;
        SP.SE.beta=NaN;
        SP.SE.kh=NaN;
        SP.SE.ke=NaN;
        
        SP.muB.mur=NaN;
        SP.muB.a=NaN;
        SP.muB.b=NaN;
        SP.muB.g=NaN;
        SP.muB.d=NaN;
        SP.muB.t=NaN;
        SP.muB.h=NaN;
        SP.muB.e=NaN;
        SP.muB.z=NaN;
        return     
        
end

% calculated mu(B) parameters
SP.muB.d = SP.muB.a./SP.muB.b;
SP.muB.t = exp(-SP.muB.b.*SP.muB.g);
SP.muB.h = SP.muB.a.*SP.muB.t;
SP.muB.e = SP.muB.t./(SP.muB.t+1);
SP.muB.z = 1./(SP.muB.t+1);

end