%% ToroidEvaluateTriangular.m
% This script evaluates a single example design, from a set of design
% specs, an example triangular excitation waveform, and design descriptive 
% variables.
% 
% Written: Vinicius C. do Nascimento
%          viniciuscn1@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% clear workspace
close all;                  % close all windows
clear;                      % clear wokspace
clc;                        % clear display

% add libs
addpath(genpath([pwd,'\lib']));

% define local variables
fn = 1;                     % figure # control variable

% define local constants
um = 1e-6;                  % convert micrometer to meter [m]
mm = 1e-3;                  % convert milimeter to meter [m]

%% Design specifications
% general specs
D.Lrqi  = 1e-3;             % minimum required incremental inductance [H]
D.kpfmx = 0.7;              % maximum allowed packing factor
D.kb    = 1.05;             % winding build factor
D.amxa  = 10;               % aspect ratio - maximum allowed
D.MLmxa = 1.5;              % inductor mass - maximum allowed [kg]
D.PLmxa = 1e2;              % inductor loss - maximum allowed [W]
D.Bmxa  = inf;              % max section spatial average flux density allowed [T]
D.nspc  = 1;                % # of parallel strands

% manufacturing limits
D.dep   = 0.8*mm;           % core protection layer thickness [m]
D.lwrmx = 9.5;              % max. wire length-toroidal winding machine [m]
D.dsmx  = inf;              % max. strand diameter (inf->no max strand diamater)

% protection layer material properties
D.ep.rho = 1.25e3;          % protection layer material density [kg/m^3]

% thermal design specifications
D.ti    = 27*um;            % assumed wire insulation thickness [m]
D.emax  = 0.01;             % thermal error criteria (K)
D.kmax  = 100;              % maximum thermal interation
D.mii   = 9;                % thermal material index of wire insulation
D.mia   = 8;                % thermal material index of air
D.micd  = 1;                % thermal material index of copper
D.hca   = 15;               % heat transfer coefficeint 
                            % core to air [W/(m^2 K)]
D.hwa   = 25;               % heat transfer coefficient from 
                            % winding bundle to air [W/(m^2 K)]
D.hsl   = 560;              % protection layer heat transfer coefficient (W/(m^2 K))  
D.gwc   = 0;                % distance from protection layer to winding bundle (m)     
D.Ta    = 25+273.15;        % ambient temperature
D.Trise = 80;               % temperature rise 
D.Twmxa = D.Ta+D.Trise;     % winding temperature limit (K)
D.Tcrmxa = D.Ta+D.Trise+50; % core temperature limit (K)

%% Coil excitation
% There are two options:
%  1) Define current waveform by setting a time vector (D.t) in [s] and 
%  a current vector (D.Ic) in [A].
%  Example:
%       fsw  = 10e3;                % switching frequency [Hz]
%       f1   = 60;                  % fundamental frequency [Hz]
%       D.t  = transpose(linspace(0,1/f1,round(fsw/f1*100)));
%       D.Ic = 5+15*cos(2*pi*f1*D.t)+1*cos(2*pi*fsw*D.t);% coil current waveform [A]
%  2) Define an array of harmonic frequencies (D.fcf) in [Hz], harmonic 
%  currents (D.Icf) in [A] and another for phases (D.phf) in [rad], such 
%  that the current can be defined in time as follows:
%  ifit = sum(D.Icf.*cos(2*pi*D.fcf*D.t'+D.phf));
%  Example:
%       D.fcf = [0,60,10e3]; % [DC,fundamental,switching]
%       D.Icf = [5,15,1];
%       D.phf = [0,0,0];

% % define current waveform applied to the coil
% fsw  = 5e3;                % switching frequency [Hz]
% f1   = 60;                  % fundamental frequency [Hz]
% % have at least 100 points for every switching cycle
% D.t  = transpose(linspace(0,1/f1,round(fsw/f1*100)));
% D.Ic = 5+15*cos(2*pi*f1*D.t)+1*cos(2*pi*fsw*D.t);% coil current waveform [A]
% D.Ic = 15*ones(length(D.t),1);

% create a waveform to evaluate losses: 10 A_dc + triangular ripple at 50%
% duty cycle. Note that the time vector only cover 1 period of the minimum
% frequency in the signal, in this case 20kHz
fsw = 20e3;     % switching frequency [Hz]
Idc = 10;       % DC current [A]
Irpp = 1;       % amplitude (0 to peak) of ripple current [A]
D.t = transpose(linspace(0,1/fsw));
D.Ic = Idc+Irpp*sawtooth(2*pi*fsw*D.t,1/2); 
% you may choose how many harmonics you want analyzed, but that increases
% the computational burden. If you do not define D.nf (# of harmonics), it
% will automatically pick the most important harmonics (greater than 1% of
% fundamental).
D.nf = 10;                  % force to evaluate 10 harmonics

% define current structure
[D,fn] = DefineCurrent(D,fn);

%% Core Losses
% This section defines how the core losses will be computed. The core
% losses calculation in this toolbox is a behaviroual model based on the 
% Steinmetz Equation framework. The parameters for calculation are defined
% in the CoreMaterialCatalog. These are the standard kh, alpha, and beta,
% with the addition of the eddy current losses term ke. Two methods to
% compute the 'hysteresis' losses are available (MSE and iGSE), and it is
% selected in D.cl.method. The frequency for core loss calculation is
% defined in D.cl.fcl in [Hz]. Notice that the flux density for either
% method needs to be computed for an evaluation waveform. The number of
% evaluation points is selected in D.cl.ncl. 
% More specific to how the core losses parameters was fitted, one can 
% select how the evaluation waveform will be created. Two options are
% available: 'sine' and 'full'. In the former, the evaluation waveform is a
% sinusoid, therefore, regardless of the waveform defined in the coil
% excitation, the core losses is calculated assuming a sine wave of
% frequency D.cl.fcl and applied current amplitude D.cl.Icl. The second
% option is 'full', where whatever waveform provided in the coil excitation
% is processed for core losses calculation. Note that, the provided coil
% excitation waveform is undesampled to D.cl.ncl points. Finally, a
% resource to speed up the evaluation is provided through D.cl.symmetry.
% This parameter allows the evaluation of just a fraction of the evaluation 
% waveform. For instance, if that is a sine wave, setting it 1/2 would 
% evaluate the losses just for the positive half of it, and account for the 
% waveform symmetry in the integral of the losses. The same is true for 
% 1/4. When symmetry is less than 1, deltaB=Bmax-Bmin is assumed to be 
% deltaB=2*Bmax. When unsure about it for D.cl.wave different than 'sine', 
% leave it set to 1. 
% If you want to neglect core losses, just set the frequency D.cl.fcl to 0.
D.cl.method   = 'MSE';      % core losses methods: 'MSE' or 'iGSE'
D.cl.fcl      = fsw;        % frequency to compute core losses [Hz]
D.cl.ncl      = 20;         % # of evaluation points for core losses
D.cl.wave     = 'sine';     % 'sine': sinusoidal approximation
                            % 'full': process full waveform
D.cl.symmetry = 1/4;        % process just a fraction of the period due to 
                            % symmetry: e.g. 1/2->positive half of sine,
                            % 1/4->positive quarter of sine
% if the D.cl.wave selected is 'sine' approximation, define amplitude of 
% the sinusoidal waveform
D.cl.Icl      = Irpp;       % amplitude of ripple current to compute core losses [Apk]

%%
% get material properties
D.mp = CoreMaterialCatalog(1,D.cl.fcl);  % store material properties in parameters struct
D.cp = conductor_catalog(1);    % conductor properties (copper)
    
%%
% simulation parameters
D.ns    = 10;               % # of sections (concentric toroidal rings)
n = 1;                      % # of permeability function multipliers:
                            % 1: constant permeability across the radius
                            % 2: affine permeability function
                            % 3: quadratic permeability function
                            % 4: cubic permeability function
                            % 5 or more: piecewise cubic hermite
                            %            permeability function
D.nfmu = ones(n,1);         % array defining order of permeability tuning
D.nobj = 2;                 % # of optimization objects
                            % 2: mass and loss
                            % 3: mass, volume and loss
D.PwInt = 'hermite';        % Arbitrary function interpolation approach:
                            % 'hermite': Piecewise Hermite Cubic Polynomial
                            % 'linear': Piecewise Linear Polynomial 

%%
% Evaluate an arbitrary design
x(1) = 150;                 % 1   # of turns
x(2) = 2.08196475719e-06;   % 2   ac*: desired conductor cross section [m^2]
x(3) = 0.01;                % 3   ra*: desired toroidal core inner air radius [m]
x(4) = 0.005;               % 4   lc: toroidal core axial length [m]
x(5) = 0.02;                % 5   dc: toroidal core radial depth [m]
x(6) = 1;                   % 6+  fmu: permeability multipliers (for permeability tuning)
% call evaluation function 
ToroidEval(x,D,fn);     