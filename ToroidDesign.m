%% ToroidDesign.m
% This script performs a toroid optimization from a set of Design specs and
% an excitation waveform defined in this file. At the end of an
% optimization, it saves a result file that can be post processed using
% ToroidPostProcess.m script.
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

% initialize local variables
fn = 1;                     % figure # control variable

% define local constants
um = 1e-6;                  % convert micrometer to meter [m]
mm = 1e-3;                  % convert milimeter to meter [m]

%% Design specifications
% general specs
D.Lrqi  = 0.75e-3;          % minimum required incremental inductance [H]
D.kpfmx = 0.7;              % maximum allowed packing factor
D.kb    = 1.05;             % winding build factor
D.amxa  = 10;               % aspect ratio - maximum allowed
D.MLmxa = 1;                % inductor mass - maximum allowed [kg]
D.PLmxa = 1e2;              % inductor loss - maximum allowed [W]
D.Bmxa  = inf;              % max section spatial average flux density allowed [T]
D.nspc  = 1;                % # of parallel strands

% manufacturing limits
D.dep   = 0.81*mm;          % epoxy layer for core protection thickness [m]
D.lwrmx = 9.5;              % max. wire length-toroidal winding machine [m]
D.dsmx  = 1.024*mm;         % max. strand diameter (AWG 18)

% epoxy protection layer material properties
D.ep.rho = 1.25e3;          % density [kg/m^3]

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
D.Ta    = 40+273.15;        % ambient temperature
D.Trise = 80+10;            % temperature rise 
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

% define the currents with your frequency domain
D.fcf = [0,1e3];    % [0 means DC, frequencies of the present harmonics]
D.Icf = [0,15];     % corresponding amplitudes [DC current, ripple current of 1 A for example]
D.phf = [0,0];      % no phase in harmonics

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
D.cl.method   = 'MSE';      % core losses methods: 'MSE' or 'iGSE'
D.cl.fcl      = 1e3;        % frequency to compute core losses [Hz]
D.cl.ncl      = 20;         % # of evaluation points for core losses
D.cl.wave     = 'sine';     % 'sine': sinusoidal approximation
                            % 'full': process full waveform
D.cl.symmetry = 1/4;        % process just a fraction of the period due to 
                            % symmetry: e.g. 1/2->positive half of sine,
                            % 1/4->positive quarter of sine
% if the D.cl.wave selected is 'sine' approximation, define amplitude of 
% the sinusoidal waveform
D.cl.Icl      = 15;       % amplitude of ripple current to compute core losses [Apk]

%%
% get material properties
D.mp = CoreMaterialCatalog(1,D.fcl);% store material properties in parameters struct
D.cp = conductor_catalog(1);        % conductor properties (copper)

%%
% simulation parameters
D.ns = 10;                  % # of sections (concentric toroidal rings)
n = 3;                      % # of permeability function multipliers:
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
% genetic algorithm parameters
ngen = 1e3;                     % # of generations
npop = 1e3;                     % size of population
obj  = 0;                       % objective to optimize
GAP=gapdefault(D.nobj,obj,npop,ngen); % setup the GA to its default config.
GAP.ev_pp = 1;                  % enable parallel loop
GAP.rp_lvl = 0;                 % disable plotting while optimize to speed it up

% gene description---------------------------------------------------------
%
%         min   max chrm  chrm      par
%         val   val type   id         #     description
GAP.gd=[1e0    1e4   3    1;  ...  %  1   # of turns
        1e-9   1e-3  3    1;  ...  %  2   ac*: desired conductor bundle area [m^2]
        1e-2   0.1   3    1;  ...  %  3   ra*: desired toroidal core inner air radius [m]
        5e-3   0.1   3    1;  ...  %  4   lc: toroidal core axial length [m]
        1e-3   0.1   3    1;  ...  %  5   dc: toroidal core radial depth [m]
        1*D.nfmu    100/D.mp.mur*D.nfmu   3*D.nfmu  1*D.nfmu]; %  6~6+D.ns   f(mu)

% perform optimization
[fP,GAS,bi,bf] = gaoptimize(@ToroidFit,GAP,D);

% save results
filename = ['Results_',num2str(D.Twmxa-273.15),'oC_',...
    'k',num2str(length(D.nfmu)),'_',...
    regexprep(regexprep(datestr(now),' ','_'),':','-')];
save(filename);