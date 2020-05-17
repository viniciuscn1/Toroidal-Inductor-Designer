%% ToroidEvaluateSine.m
% This script evaluates a single example design, from a set of design
% specs, an example sinusoidal excitation waveform, and design descriptive 
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

% initialize local variables
fn = 1;                     % figure # control variable

% define local constants
um = 1e-6;                  % convert micrometer to meter [m]
mm = 1e-3;                  % convert milimeter to meter [m]

%% Design specifications
% general specs
D.Lrqi  = 1e-3;             % minimum required incremental inductance [H]
D.kpfmx = 0.7;              % maximum allowed packing factor
D.kb    = 1.05;             % winding build factor
D.amxa  = 50;               % aspect ratio - maximum allowed
D.MLmxa = 1.5;              % inductor mass - maximum allowed [kg]
D.PLmxa = 1e2;              % inductor loss - maximum allowed [W]
D.Bmxa  = inf;              % max section spatial average flux density allowed [T]
D.nspc  = 1;                % # of parallel strands

% manufacturing limits
D.dep   = 0.8*mm;           % epoxy layer for core protection thickness [m]
D.lwrmx = 9.5;              % max. wire length-toroidal winding machine [m]
D.dsmx  = inf;              % max. strand diameter (inf->no max diameter)

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

% define the currents with your frequency domain
D.fcf = [0,1e3];  % [DC,switching frequency]
D.Icf = [0,15]; % [DC current, ripple current of 1 A for example]
D.phf = [0,0];      % no phase in harmonics

% current to compute core losses
D.Icl= 15;         % amplitude of ripple current to compute core losses [Apk]
D.fcl= 1e3;        % frequency to compute core losses [Hz]

% define current structure
[D,fn] = DefineCurrent(D,fn);

%%
% get material properties
D.mp = CoreMaterialCatalog(1,D.fcl);% store material properties in parameters struct
D.cp = conductor_catalog(1);        % conductor properties (copper)

%%
% simulation parameters
D.ns = 10;                  % # of sections (concentric toroidal rings)
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