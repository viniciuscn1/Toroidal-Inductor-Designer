function [f] = ToroidFit(x,D,fn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitness function to evaluate the toroidal core inductor metrics         %
%                                                                         %
% Written by Vinicius Nascimento                                          %
% viniciuscn1@gmail.com
%                                                                         %
% Current revision date: 01-10-2020                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Call: [f] = ToroidFit(x,D,fn);
%       [f] = ToroidFit(x,D);
% 
% Inputs: 
% 
%  x        = optimization variables (genes)
%   x(1)    = desired # of turns
%   x(2)    = ac*: desired conductor cross section [m^2]
%   x(3)    = ra*: desired toroidal core inner air radius [m]
%   x(4)    = lc: toroidal core axial length [m]
%   x(5)    = dc: toroidal core radial depth [m]
%   x(6:end)= fmu: permeability multipliers (for permeability tuning)
%  D        = design parameter structure
%  fn       = figure number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% determine mode: evaluation or fitness
if (nargin>2)       % evaluation mode
    ev = true;      % set evaluation flag on
else                % fitness mode
    ev = false;     % set evaluation flag off
end

%% Optimization variables and local constants and variables definition
% define local constants
NC = 10;                        % total number of constraints
c  = zeros(NC,1);               % constraints flags array

% setup parameters variable
P.mp = D.mp;                    % initialize core material properties
P.cp = D.cp;                    % initialize conductor material properties

% extract optimization variables
P.N = round(x(1));              % # of turns 
WP = AwgRoundMax(sqrt(x(2)/(pi*D.nspc)),D.dsmx);  % get wire strand gauge properties
P.as = WP.area;                 % conductor strand area [m^2]
P.ac = P.as*D.nspc;             % conductor area [m^2]
P.desc = WP.desc;               % AWG description
P.ra = x(3);                    % inner window air radius [m]
P.lc = x(4);                    % toroidal core axial length [m]
P.dc = x(5);                    % toroidal core radial depth (ro-ri) [m]
fmu = x(6:end);                 % permeability multiplier datapoints
P.fmu = fmu(:);                 % enforce column vector

if length(P.fmu)==length(D.nfmu)% check if the # of datapoints is ok
    if length(P.fmu)==1         % constant permeability throughout toroid
        P.nmu = 1;              % number of permeability sections
    else
        P.nmu = length(P.fmu)-1;% number of permeability sections
    end
else
    error('Number of permeability multipliers is not compatible to the number of genes.');
end

%% Geometry
% calculate dependent geometry (remaining dimensions)

% get actual conductor radius from the selected in the AWG catalog
P.rs = sqrt(P.as/pi);           % get conductor strand radius [m]
P.rc = sqrt(P.ac/pi);           % get conductor bundle equivalent radius [m]
% compute # of coil layers based on the # of turns (inner side of coil)
P.Nli = 1;	% # of coil layers based on the # of turns (inner side of coil)
Ncap(P.Nli) = pi/(asin(P.rs*D.kb/(P.ra+(2*P.Nli-1)*P.rs*D.kb)));
while sum(Ncap)<P.N*D.nspc
    P.Nli = P.Nli+1;
    Ncap(P.Nli) = pi/(asin(P.rs*D.kb/(P.ra+(2*P.Nli-1)*P.rs*D.kb)));
end

% for the inner coil, the number of layers was determined by adding layers
% from inside out. However, the actual inductor has the conductors filling
% the for the outermost layer of the inner coil inwards. Since the
% outermost layer has a greater conductor hosting capacity, the innermost
% lays might be empty, so this has to be corrected as follows
if sum(cumsum(fliplr(Ncap))>P.N*D.nspc)>1
    P.Nli = find(cumsum(fliplr(Ncap))>P.N*D.nspc,1,'first');
    P.ra  = P.ra+2*P.rs*D.kb*(length(Ncap)-P.Nli);
end

% compute radial depth of inner winding bundle dwi
P.dwi = 2*P.rs*D.kb*P.Nli;

% compute toroidal core dimensions 
P.rei = P.ra+P.dwi;                 % radius to inner side of epoxy [m]
P.rci = P.rei+D.dep;                % core inner radius [m]
P.rco = P.rci+P.dc;                 % core outer radius [m]
P.rwo = P.rco+D.dep;                % radius to outer winding bundle [m]

% compute # of coil layers based on the # of turns (outer side of coil)
P.Nlo = 1;	% # of coil layers based on the # of turns (outer side of coil)
Ncap = pi/(asin(P.rs*D.kb/(P.rco+(2*P.Nlo-1)*P.rs*D.kb)));
while Ncap<P.N*D.nspc
    P.Nlo = P.Nlo+1;
    Ncap = Ncap + pi/(asin(P.rs*D.kb/(P.rwo+(2*P.Nlo-1)*P.rs*D.kb)));
end

% compute radial depth of outer winding bundle dwo
P.dwo = 2*P.rs*D.kb*P.Nlo;

% compute remaining geometric features
P.ro = P.rwo+P.dwo;                 % toroidal inductor outer radius [m]

% discretization or sectionalization support variables
P.rsec = transpose(linspace(P.rci,P.rco,D.ns+1));% radii breakpoints for each concentric section
P.rsecm= (P.rsec(1:end-1)+P.rsec(2:end))/2;% mean radii of each section [m]
P.dsec = diff(P.rsec);              % radial depth of each section [m]
P.Acsec= P.dsec*P.lc;               % cross sectional area of each section [m^2]
P.Assec= pi*P.dsec.*(P.rsec(1:D.ns)+P.rsec(2:D.ns+1)); % surface area of each section [m^2]
P.Vsec = P.Assec*P.lc;              % volume of each section [m^3]

%%
% Permeability function

% first create the permeability function sectionalization
P.rmu = transpose(linspace(P.rci,P.rco,P.nmu+1));

% permeability function variables  
if length(P.fmu)==1         % constant permeability across radius 
    P.mur = mkpp(P.rmu,[0,P.fmu]);
elseif length(P.fmu)==2     % affine permeability across radius
    P.mur = mkpp(P.rmu,[diff(P.fmu)/diff(P.rmu),P.fmu(1)]);
elseif isempty(P.fmu)       % no permeability multiplier defined (error)
    error('Permeability function %s is invalid\n.');    
else                        % arbitrary permeability function
    if isfield(D,'PwInt')
        if strcmp(D.PwInt,'hermite')
            % fit a piecewise cubic function to the datapoints provided
            P.mur = pchip(P.rmu,P.fmu);
        elseif strcmp(D.PwInt,'linear')
            % fir a piecewise linear function to the datapoints provided
            P.mur = mkpp(P.rmu,[diff(P.fmu)./diff(P.rmu),P.fmu(1:end-1)]);
        else
            error(['The piecewise interpolation apprpach %s is invalid. ',...
            'The options are ''Hermite'' and ''Linear''.'],D.PwInt);
        end
    else
        error(['Please define the field PwInt to the D structure. ',...
            'The options are ''Hermite'' and ''Linear''.']);
    end
end

%% Volume, mass and packing factor
% compute coil area, volume and resistance for each section
% herein, the coil is divided in the following sections:
% - inner side;
% - outer side;
% - top and bottom sections;
% - inner corner (2x transitions from inner side to top or bottom)
% - outer corner (2x transitions from outer side to top or bottom)

% compute coil cross sectional area for each section ----------------------

% compute inner side coil cross sectional area 
% Acli = pi*(rei^2-ra^2) = pi*(rei-ra)*(rei+ra) = pi*dwi*(rei+ra)
Acli = pi*P.dwi*(P.rei+P.ra);
% compute outer side coil cross sectional area
% Aclo = pi*(ro^2-rwo^2)=pi*dwo*(ro+rwo)
Aclo = pi*P.dwo*(P.ro+P.rwo);
% compute average area of top (upper) or bottom (lower) sides (horizontal
% section)
% first, the height of the bundle above ri and ro have to be determined. 
% The width is defined by the circumference around inner and outer radius
% of the toroidal core
Nwi = floor(2*pi*P.rci/(2*P.rs*D.kb));  % # of conductors in row above ri
Nwo = floor(2*pi*P.rco/(2*P.rs*D.kb));  % # of conductors in row above ro    
% compute # of vertical conductor layers above ri and ro
Ndi = ceil(D.nspc*P.N/Nwi);             % # of conductors in column above ri
Ndo = ceil(D.nspc*P.N/Nwo);             % # of conductors in column above ro
% compute vertical height above ri and ro
P.hri = Ndi*2*P.rs*D.kb;                % conductor height above ri [m]
P.hro = Ndo*2*P.rs*D.kb;                % conductor height above ro [m]
% divide the coil in ns sections as well (similar to the core)
% inner and outer height of each coil section [m]
P.hws = transpose(linspace(P.hri,P.hro,D.ns+1));
P.hwm = (P.hws(1:end-1)+P.hws(2:end))/2;
% note that Aclh(r) = 2*pi*r*h, where h(ri) = hri and h(ro) = hro
% Aclh(r) = 2*pi*r/(ro-ri)*[(hro-hri)*r + hri*ro - hro*ri] 
% solving the incremental resistance problem dR=dr/(sigma.A(r))
% Aclh = 2*pi*(W.hri*P.rco-W.hro*P.rci)/log(P.rco*W.hri/(P.rci*W.hro));
% compute the coil equivalent area of horizontal region for each section
% equally divided according to the core sectionalization
Acls = (2*pi*(P.hws(1:D.ns).*P.rsec((1:D.ns)+1)-...
        P.hws((1:D.ns)+1).*P.rsec(1:D.ns))./...
        log(P.rsec((1:D.ns)+1).*P.hws(1:D.ns)./...
        (P.rsec(1:D.ns).*P.hws((1:D.ns)+1))));
% compute transition section average cross sectional areas
% in order to account for the variation of the bundle thickness around the
% corner l(theta), consider taking the average between dwi and hri
Aclic = pi*(P.dwi+P.hri)/2*sqrt((2*P.rci)^2-(2*D.dep+(P.dwi+P.hri)/2)^2);
Acloc = pi*(P.dwo+P.hro)/2*sqrt((2*P.rco)^2-(2*D.dep+(P.dwo+P.hro)/2)^2);

% compute average packing factor for each section -------------------------
Acd  = D.nspc*P.N*P.as;     % conductor total cross sectional area [m^2]
kpfi = Acd/Acli;            % inner side packing factor 
kpfo = Acd/Aclo;            % outer side packing factor 
kpfic= Acd/Aclic;           % inner corners average packing factor
kpfoc= Acd/Acloc;           % outer corners average packing factor
kpfs = Acd./Acls;           % packing factor of each section of the winding 
                            % horizontal region
% merge all sections packing factors
kpf = [kpfi;kpfo;kpfs;kpfic;kpfoc];

% compute coil volume for each coil section -------------------------------
% compute the volume of each horizontal section
P.Vcls = 2*pi./(P.rsec(2:end)-P.rsec(1:end-1)).*((P.hws(2:end)-P.hws(1:end-1))/3.*...
    (P.rsec(2:end).^3-P.rsec(1:end-1).^3)+(P.rsec(2:end).^2-P.rsec(1:end-1).^2)/2.*...
    (P.rsec(2:end).*P.hws(1:end-1)-P.hws(2:end).*P.rsec(1:end-1)));
P.Vcli = Acli*P.lc;         % coil inner side volume [m^3]
P.Vclo = Aclo*P.lc;         % coil outer side volume [m^3]
Vclic= pi^2*P.rci*(D.dep*(P.hri+P.dwi)+P.dwi*P.hri)/2;% inner corner vol. [m^3]
Vcloc= pi^2*P.rco*(D.dep*(P.hro+P.dwo)+P.dwo*P.hro)/2;% outer corner vol. [m^3]

% compute total coil volume [m^3]
P.Vcl = P.Vcli+P.Vclo+2*sum(P.Vcls)+2*Vclic+2*Vcloc;

% compute conductor volume of each coil section ---------------------------
P.Vcdi = kpfi*P.Vcli;           % conductor inner side volume [m^3]
P.Vcdo = kpfo*P.Vclo;           % conductor outer side volume [m^3]
P.Vcdic= kpfic*Vclic;           % conductor inner corners volume [m^3]
P.Vcdoc= kpfoc*Vcloc;           % conductor outer corners volume [m^3]
P.Vcds = kpfs.*P.Vcls;          % conductor sections [m^3]

% compute total conductor volume [m^3]
P.Vcd = P.Vcdi+P.Vcdo+2*sum(P.Vcds)+2*P.Vcdic+2*P.Vcdoc;

% compute total wire length [m]
P.lwr = P.Vcd/P.ac;

% compute coil resistance -------------------------------------------------
P.Ri = P.Vcli*P.N^2/(kpfi*Acli^2*P.cp.sigma0);
P.Ro = P.Vclo*P.N^2/(kpfo*Aclo^2*P.cp.sigma0);
P.Ric= Vclic*P.N^2/(kpfic*Aclic^2*P.cp.sigma0);
P.Roc= Vcloc*P.N^2/(kpfoc*Acloc^2*P.cp.sigma0);
P.Rh = P.Vcls.*P.N^2./(kpfs.*Acls.^2*P.cp.sigma0);
% compute total coil series resistance [Ohm]
P.Rs = P.Ri+P.Ro+2*sum(P.Rh)+2*P.Ric+2*P.Roc;

%%
% compute total mass
P.Vep = pi*P.dc*(P.rco+P.rci)*D.dep*2 + ... epoxy top and bottom 
        pi*D.dep*(P.rei+P.rci)*P.lc +   ... epoxy inner side
        pi*D.dep*(P.rco+P.rwo)*P.lc +   ... epoxy outer side
        (pi*D.dep)^2*(P.rci+P.rco);     ... epoxy corners
P.Mep = P.Vep*D.ep.rho;             % mass of epoxy protection layer [kg]
% compute toroidal core volume [m^3]
P.Vcr = pi*P.dc*P.lc*(P.rco+P.rci); % toroidal core volume [m^3]
P.Mcr = P.Vcr*P.mp.rho;             % toroidal core mass [kg]
P.Mcd = P.Vcd*P.cp.row;             % total conductor mass [kg]
P.ML  = P.Mcd+P.Mcr+P.Mep;          % total inductor mass [kg]

% compute peak current density
P.Jpk = D.Ipk/P.ac;                 % conductor peak AC current density [A/m^2]

%%
% compute circumscribing volume and aspect ratio
P.hL = P.lc+2*P.hri+2*D.dep;            % total height (axial length) [m]
P.dL = 2*P.ro;                          % total diameter [m]
P.aL = max([P.hL,P.dL])/min([P.hL,P.dL]);% aspect ratio
P.VL = pi*(P.dL/2)^2*P.hL;              % total circumscribing volume [m^3]

%%
% check constraints
c(1) = lte(P.lwr,D.lwrmx);          % constrain coil max wire length
c(2) = lte(max(kpf),D.kpfmx);       % constrain max packing factor
c(3) = lte(P.aL,D.amxa);            % constrain max aspect ratio
c(4) = lte(P.ML,D.MLmxa);           % constrain max inductor mass
% test imposed constraints 
CS = sum(c);                        % constraints satisfied
CI = 4;                             % constraints imposed
if (CS<CI)&&(~ev)                   % if any constraint was not satisfied
   f=eps*ones(D.nobj,1)*(CS-NC)/NC; % return a negative fitness proportional
   return;                          % to the number of satisfied constraints
end

%%
% compute device's incremental inductance [H]
P    = ComputeLeakage(D,P); % compute leakage inductance [H] (neglect corners)
r    = transpose(1./linspace(1/P.rci,1/P.rco,D.ns));% 1/r spacing
Ip   = [0.99,1]*D.Ipk;      % incremental current at peak value [A] (row vector)
B    = SolveB(P,Ip,r);      % solve for flux density [T]
lamc = P.N*trapz(r,B)*P.lc; % core flux linkage [Vs]
P.Lm = diff(lamc)./diff(Ip);% core magnetizing inductance at full load [H]
P.Linc = P.Lm+P.Llk;        % inductor incremental inductance at full load [H]
c(5) = gte(P.Linc,D.Lrqi);  % enforce min. L requirement is met
c(6) = lte(max(B(:,end)),D.Bmxa);   % limit maximum allowed flux density
% test imposed constraints 
CS = sum(c);                        % constraints satisfied
CI = 6;                             % constraints imposed
if (CS<CI)&&(~ev)                   % if any constraint was not satisfied
   f=eps*ones(D.nobj,1)*(CS-NC)/NC; % return a negative fitness proportional
   return;                          % to the number of satisfied constraints
end

%% 
% compute core losses
P = ComputeCoreLosses(P,D);

%% Heat Transfer Analysis
% Thermal Equivalent Circuit (TEC)
% solve thermal equivalent circuit: computes mean and peak temperatures and
% updates winding losses according to the temperature
[T,P.Pw,P.Pp,c(7),Rdc] = ToroidTEC(D,P,P.Plcsec);

% add all losses
P.Pt = P.Plc+P.Pw+P.Pp;

% check constraints on losses and maximum temperature
c(8) = lte(P.Pt,D.PLmxa);     % constraint on maximum losses against limit
c(9) = lte(T.Twmx,D.Twmxa);   % check maximum winding temperature vs limit
c(10)= lte(T.Tcrmx,D.Tcrmxa); % check maximum core temperature vs limit
% test constraints 
CS = sum(c);                        % constraints satisfied
CI = length(c);                     % constraints imposed
if (CS<CI)&&(~ev)                   % if any constraint was not satisfied
   f=eps*ones(D.nobj,1)*(CS-NC)/NC; % return a negative fitness proportional
   return;                          % to the number of satisfied constraints
end % end of constraints test

%%
% compute fitness
if D.nobj == 3                      % # of objectives selected = 3
    f = [1/P.ML; 1/P.VL; 1/P.Pt];   % include volume minimization 
elseif D.nobj == 2                  % # of objectives selected = 2
    f = [1/P.ML; 1/P.Pt];           % only optimizes for mass and losses
else                                % other # of optimization objectives is invalid
    error('The number of objectives %d is invalid.\n',D.nobj);
end

%%
if (nargin>2)
    np = 1e2;                       % # of points for plotting
    r = linspace(P.rci,P.rco,np);   % core radii for plotting
    perf = ppval(P.mur,r);          % fitted permeability function
    if (fn>0)
        %%
        % display optimization results for selected point in pareto
        % core
        fprintf('Toroidal Core Data: --------------------\n');
        fprintf('Core material: %s \n',P.mp.desc);
        fprintf('Core mass: %.3f [kg]\n',P.Mcr);
        fprintf('Finished inner radius: %.3f [mm]\n',P.ra*1e3);
        fprintf('Inner core radius: %.3f [mm]\n',P.rci*1e3);
        fprintf('Outer core radius: %.3f [mm]\n',P.rco*1e3);
        fprintf('Outer winding radius: %.3f [mm]\n',(P.rco+P.dwo)*1e3);
        fprintf('Axial length (height): %.3f [mm]\n',P.lc*1e3);
        fprintf('Total core loss: %.3f [W]\n',P.Plc);
        fprintf('Maximum core temperature vs constraint: %.3f / %.3f [oC]\n',...
            T.Tcrmx-273.15,D.Tcrmxa-273.15);
        % winding
        fprintf('\nToroidal Winding Data: --------------------\n');
        fprintf('Material: %s\n',P.cp.desc);
        fprintf('AWG: %s\n',P.desc);
        fprintf('Total wire length vs constraint: %.3f / %.3f [m]\n',...
            P.lwr,D.lwrmx);
        fprintf('Coil mass: %.3f [kg]\n',P.Mcd);
        fprintf('Conductor area: %.3f [mm^2]\n',P.ac*1e6);
        fprintf('Conductor strand radius: %.3f [mm]\n',P.rs*1e3);
        fprintf('Conductor bundle equivalent radius: %.3f [mm]\n',P.rc*1e3);
        fprintf('Number of turns: %d \n',P.N);
        fprintf('Number of winding layers (inner side): %d \n',P.Nli);
        fprintf('Number of winding layers (outer side): %d \n',P.Nlo);
        fprintf('Total ohmic loss: %.3f [W]\n',P.Pw);
        fprintf('Total proximity losses: %.3f [W]\n',P.Pp);
        fprintf('Total winding losses: %.3f [W]\n',P.Pp+P.Pw);
        fprintf('Maximum winding temperature vs constraint: %.3f / %.3f [oC]\n',...
            T.Twmx-273.15,D.Twmxa-273.15);
        % optimization metrics and specifications
        fprintf('\nMetrics: --------------------\n');
        fprintf('Total mass vs allowed: %.3f / %.3f [kg]\n',P.ML,D.MLmxa);
        fprintf('Total volume: %.3f [L]\n',P.VL*1e3);
        fprintf('Total loss vs allowed: %.3f / %.3f [W]\n',P.Pt,D.PLmxa);
        fprintf('Magnetizing Incremental Inductance : %.3f [mH]\n',...
            P.Lm*1e3);
        fprintf('Leakage Inductance : %.3f [mH]\n',...
            P.Llk*1e3);
        fprintf('Total Incremental Inductance vs required: %.3f / %.3f [mH]\n',...
            P.Linc*1e3,D.Lrqi*1e3);
        fprintf('Resistance at Ambient Temperature: %.3f [mOhm]\n',P.Rs*1e3);
        fprintf('Resistance at Operating Temperature: %.3f [mOhm]\n',Rdc*1e3);
        %%
        % Plotting
        
        % plot permeability multipliers and fitted permeability function
        set(figure(fn),'color','white'); fn=fn+1;
        plot(r,perf,'-',P.rmu,P.fmu,'o');
        grid on; grid minor;
        xlabel('$r$, m','interpreter','latex');
        ylabel('$f_{\mu}(r)$','interpreter','latex');
        title('Permeability multiplier function','interpreter','latex');
        legend({'Permeability multiplier function',...
                'Permeability multiplier datapoints'},...
                'interpreter','latex','location','best');
        ylim([0,max(P.fmu)*1.2]);
        
        % compute lambda-i curve
        set(figure(fn),'color','white'); fn=fn+1;
        hold on;
        Ip = transpose(linspace(0,1.2*D.Ipk,30)); % range of applied current
        [~,idx] = min(abs(Ip-D.Ipk));   % find index of current sweep closest to peak
        Bsp     = SolveB(P,Ip,r);       % B for current sweep [T]
        Bpk     = SolveB(P,D.Ipk,r);    % B for peak current [T]
        lamsp   = P.N*trapz(r(:),Bsp)*P.lc+P.Llk*Ip(:)'; % sweep current flux linkage [Vs]
        lampk   = P.N*trapz(r(:),Bpk)*P.lc+P.Llk*D.Ipk;  % peak current flux linkage [Vs]
        mupk    = fmuB(P.mp,Bpk,perf(:));
        
        % plot lambda-i curve
        plot(Ip,lamsp,'o-',D.Ipk,lampk,'x-','LineWidth',2);
        hold on;
        patch([Ip(1:idx-1);D.Ipk;0],[lamsp(1:idx-1)';lampk;lampk],...
            [0.5,0.5,0.5],'FaceAlpha',0.5,'EdgeAlpha',0);
        cap = {'Analytical','Peak Current','Field Energy Stored'};
        hold off;
        grid on; grid minor;
        Wf = trapz(lamsp(1:idx),Ip(1:idx));   % compute field energy [J]
        title(['$\lambda-i$ curve - Field Energy $W_f$ = ',...
            num2str(Wf),' J'],'interpreter','latex');
        xlabel('$I_p$, A','interpreter','latex');
        ylabel('$\lambda_c$, Vs','interpreter','latex');
        legend(cap,'interpreter','latex','location','best');
        
        % plot mur(r,B)
        Bv = transpose(linspace(0,max(Bpk)*1.2,np)); % range B from 0 to Bmx
        [rr,bb] = meshgrid(r,Bv);           % grid of radius and B values
        
        % compute relative permeability
        kmu  = ppval(P.mur,r);              % tuning parameter
        mura = zeros(size(rr));             % magnetic permeability of 3D plot
        for k = 1:np
            mura(:,k) = fmuB(P.mp,Bv(k),kmu)/(4e-7*pi);
        end
        % plot affine relative permeability mur(r,B) - 3D
        set(figure(fn),'color','white'); fn=fn+1;
        surf(rr,bb,mura','EdgeColor','None');
        hold on;
        plot3(r,Bpk,mupk/(4e-7*pi),'-+','LineWidth',2,'color',[0.4940 0.1840 0.5560]);
        xlabel('$r$, m','interpreter','latex');
        ylabel('$B$, T','interpreter','latex');
        zlabel('$\mu_r(r,B)$','interpreter','latex');
        title({['Relative permeability as a function of radius and flux density $\mu_r(r,B)$'],...
            },'interpreter','latex');
        cb = colorbar('location','eastoutside');
        grid on; grid minor;
        view(135,45);
        cblim = get(cb,'Limits');
        set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
        cb.Ruler.TickLabelFormat = '%.2f';
        set(gca,'Xdir','reverse');
        
        % plot 3D cross sectional view
        fn = PlotCrossSection(fn,P,D);
        % plot 3D cross sectional view with relative magnetic permeability
        fn = PlotCrossSection(fn,P,D,Bpk);
        
        % plot temperature profile if TEC converged
        if c(7)==1
            fn = PlotTemperature(fn,P,D,T);
        end
        
        % compute analytical solution for flux density and permeability
        
        % plot 2D permeability profile
        B2D = transpose(repmat(Bpk,[1,np]));
        mur2D = transpose(repmat(mupk/(4e-7*pi),[1,np]));
        q = linspace(0,2*pi,np);
        [rr,qq] = meshgrid(r,q);
        set(figure(fn),'color','white'); fn=fn+1;
        surf(rr.*cos(qq),rr.*sin(qq),mur2D,'EdgeColor','None');
        xlabel('$x$, [m]','interpreter','latex');
        ylabel('$y$, [m]','interpreter','latex');
        title('Relative permeability $\mu_r$ across toroidal core',...
            'interpreter','latex');
        view(0,90);
        cb = colorbar('location','eastoutside');
        cblim = get(cb,'Limits');
        set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
        cb.Ruler.TickLabelFormat = '%.2f';
        hold on;
        daspect([1 1 1]);
        axis([-1,1,-1,1].*1.2*P.ro);
        PlotConductors(P,P,D);
        % plot epoxy protection layers
        PlotCircle(0,0,P.rei,P.rci,([3, 86, 31]./255),100);
        PlotCircle(0,0,P.rco,P.rwo,([3, 86, 31]./255),100);
        
        % plot 2D flux density profile
        set(figure(fn),'color','white'); fn=fn+1;
        plt = pcolor(rr.*cos(qq),rr.*sin(qq),B2D);
        set(plt,'FaceColor','interp','EdgeColor','interp');
        xlabel('$x$, [m]','interpreter','latex');
        ylabel('$y$, [m]','interpreter','latex');
        title('Flux density $B$ [T] across toroidal core',...
            'interpreter','latex');
        view(0,90);
        cb = colorbar('location','eastoutside');
        cblim = get(cb,'Limits');
        cblim(1) = 0;
        caxis([cblim(1),cblim(2)]);
        set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
        cb.Ruler.TickLabelFormat = '%.2f';
        title(cb,'$B$ [T]','interpreter','latex');
        hold on;
        daspect([1 1 1]);
        axis([-1,1,-1,1].*1.2*(P.rco+P.dwo));
        PlotConductors(P,P,D);
        % plot epoxy protection layers
        PlotCircle(0,0,P.rei,P.rci,([3, 86, 31]./255),100);
        PlotCircle(0,0,P.rco,P.rwo,([3, 86, 31]./255),100);
        hold on;
        [Bpkmx,idxmx] = max(Bpk);
        txtmx = [' ',num2str(Bpkmx,3)];
        text(r(idxmx)*cos(0),r(idxmx)*sin(0),txtmx);
        pmx = plot(r(idxmx)*cos(0),r(idxmx)*sin(0),'rx','LineWidth',2);
        [Bpkmn,idxmn] = min(Bpk);
        txtmn = [' ',num2str(Bpkmn,3)];
        text(r(idxmn)*cos(pi),r(idxmn)*sin(pi),txtmn);
        pmn = plot(r(idxmn)*cos(pi),r(idxmn)*sin(pi),'bx','LineWidth',2);
        legend([pmx,pmn],{'max','min'},'interpreter','latex');
        
        % plot 3D
        % affine permeability function
        [xi,yi,zi] = cylinder(P.rci,np);     % inner side
        [xo,yo,zo] = cylinder(P.rco,np);     % outer side
        set(figure(fn),'color','white'); fn=fn+1;
        surf(xi,yi,P.lc*zi,mur2D(1,1)*ones(size(zi)),'EdgeColor','None');
        hold on;
        surf(xo,yo,P.lc*zo,mur2D(end,end)*ones(size(zo)),'EdgeColor','None');
        surf(rr.*cos(qq),rr.*sin(qq),P.lc*ones(size(rr)),mur2D,'EdgeColor','None');
        surf(rr.*cos(qq),rr.*sin(qq),zeros(size(rr)),mur2D,'EdgeColor','None');
        hold off;
        daspect([1 1 1]);
        grid on; grid minor;
        xlabel('$x$, [m]','interpreter','latex');
        ylabel('$y$, [m]','interpreter','latex');
        zlabel('$z$, [m]','interpreter','latex');
        title({ ['Relative permeability $\mu_r$ across toroidal core'],...
                ['3D plot']},'interpreter','latex');
        camlight;
        cb = colorbar('location','eastoutside');
        cblim = get(cb,'Limits');
        set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
        cb.Ruler.TickLabelFormat = '%.2f';
        title(cb,'$\mu_r$','interpreter','latex');
        colormap;
        shading interp;
        lighting gouraud;

        % flux density
        set(figure(fn),'color','white'); fn=fn+1;
        surf(xi,yi,P.lc*zi,B2D(1,1)*ones(size(zi)),'EdgeColor','None');
        hold on;
        surf(xo,yo,P.lc*zo,B2D(end,end)*ones(size(zo)),'EdgeColor','None');
        surf(rr.*cos(qq),rr.*sin(qq),P.lc*ones(size(rr)),B2D,'EdgeColor','None');
        surf(rr.*cos(qq),rr.*sin(qq),zeros(size(rr)),B2D,'EdgeColor','None');
        hold off;
        daspect([1 1 1]);
        grid on; grid minor;
        xlabel('$x$, [m]','interpreter','latex');
        ylabel('$y$, [m]','interpreter','latex');
        zlabel('$z$, [m]','interpreter','latex');
        title({['Flux density $B$ across toroidal core'],...
            ['3D plot']},'interpreter','latex');
        camlight;
        cb = colorbar('location','eastoutside');
        cblim = get(cb,'Limits');
        cblim(1) = 0;
        caxis([cblim(1),cblim(2)]);
        set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
        cb.Ruler.TickLabelFormat = '%.2f';
        title(cb,'$B$ [T]','interpreter','latex');
        colormap;
        shading interp;
        lighting gouraud;
    end
    clear f;
    f.c = c;             % constraints
    f.ML= P.ML;          % mass
    f.PL= P.Pt;          % power losses
    f.Pw= P.Pw;          % ohmic and skin effect losses
    f.Pp= P.Pp;          % proximity losses
    f.Pc= P.Plc;           % core losses
    f.J = P.Jpk;         % peak current density
    f.par = P;           % toroid parameters
    f.W = P;             % winding information
    f.hL= P.hL;            % height
    f.dL= P.dL;            % diameter
    f.VL= P.VL;            % volume
    f.aL = P.aL;           % aspect ratio
    f.Twmx = T.Twmx;     % max winding temperature
    f.perf = perf;       % permeability function
    f.rperf= r;          % radius
    f.lwr  = P.lwr;      % wire length
end % end if (nargin>2)
    
end
%%
function h = PlotConductors(T,W,D)
hold on;                    % use current figure
%%
% plot inner layers of conductors
layer = W.Nli;              % initialize layer counter
kcl = 1;                    % initialize conductor per layer counter
% compute conductor capacity at current layer
Ncap = floor(pi/(asin(W.rs*D.kb/(T.ra+(2*layer-1)*W.rs*D.kb))));
if layer==1
    qcl = 2*pi/(W.N*D.nspc);         % compute angle between conductor centers
else
    qcl = 2*pi/Ncap;        % compute angle between conductor centers
end
for kc = 1:W.N*D.nspc             % loop though all conductors
    % compute conductor center x and y coordinates
    xc = (T.ra+W.rs*D.kb*(2*layer-1))*cos(qcl*(kcl-1)+qcl/2*mod(layer,2));
    yc = (T.ra+W.rs*D.kb*(2*layer-1))*sin(qcl*(kcl-1)+qcl/2*mod(layer,2));
    % plot conductor
    PlotCircle(xc,yc,0,W.rs,'r'); 
    PlotCircle(xc,yc,W.rs*0.8,W.rs,'k'); 
    % check counter
    if (kcl==Ncap)  % if conductor per layer count reaches layer capacity
        kcl = 1;            % reset conductor per layer count
        layer = layer-1;    % increment layer counter
        % compute new layer conductor capacity
        Ncap = floor(pi/(asin(W.rs*D.kb/(T.ra+(2*layer-1)*W.rs*D.kb))));
        if layer==1
            qcl = 2*pi/(W.N*D.nspc-kc);% compute angle between conductor centers
        else
            qcl = 2*pi/Ncap;    % compute angle between conductor centers
        end
    else            % if layer hasn't been filled of conductors yet
        kcl = kcl+1;        % increment conductor per layer counter
    end
end
%%
% plot outer layers of conductors
layer = 1;                  % initialize layer counter
kcl = 1;                    % initialize conductor per layer counter
% compute conductor capacity at current layer
Ncap = floor(pi/(asin(W.rs*D.kb/(T.rwo+(2*layer-1)*W.rs*D.kb))));
if W.Nlo==1
    qcl = 2*pi/(W.N*D.nspc);         % compute angle between conductor centers
else
    qcl = 2*pi/Ncap;        % compute angle between conductor centers
end
for kc = 1:W.N*D.nspc              % loop though all conductors
    % compute conductor center x and y coordinates
    xc = (T.rwo+W.rs*D.kb*(2*layer-1))*cos(qcl*(kcl-1)+qcl/2*mod(layer,2));
    yc = (T.rwo+W.rs*D.kb*(2*layer-1))*sin(qcl*(kcl-1)+qcl/2*mod(layer,2));
    % plot conductor
    PlotCircle(xc,yc,0,W.rs,'r'); 
    PlotCircle(xc,yc,W.rs*0.8,W.rs,'k'); 
    % check counter
    if (kcl==Ncap)  % if conductor per layer count reaches layer capacity
        kcl = 1;            % reset conductor per layer count
        layer = layer+1;    % increment layer counter
        % compute new layer conductor capacity
        Ncap = floor(pi/(asin(W.rs*D.kb/(T.rwo+(2*layer-1)*W.rs*D.kb))));
        if layer==W.Nlo
            qcl = 2*pi/(W.N*D.nspc-kc);% compute angle between conductor centers
        else
            qcl = 2*pi/Ncap;    % compute angle between conductor centers
        end
    else            % if layer hasn't been filled of conductors yet
        kcl = kcl+1;        % increment conductor per layer counter
    end
end
hold off;                   % release figure
end
%%
function h = PlotCircle(xc,yc,ri,ro,color,np)
hold on;                        % use current figure
if nargin<6
    np = 50;                    % # of points
end
q = linspace(0,2*pi,np);        % theta ranging from 0 to 2*pi [rad]
x = [ro*cos(q),ri*cos(q)]+xc;   % x coordinates [m]
y = [ro*sin(q),ri*sin(q)]+yc;   % y coordinates [m]
h = fill(x,y,color,'EdgeColor','none');% plot conductor
hold off;                       % release figure
end
%%
function fn = PlotCrossSection(fn,P,D,Bpk)

% define local constants
np = 250;

set(figure(fn),'color','white'); fn=fn+1;
% plot cross section
yc = P.lc/2;
y1q = linspace(yc,yc+D.dep+P.hro,np);% y coord.(1st) [m]
x1q = P.rco+(D.dep+P.dwo)*sqrt(1-min(((y1q-yc)./(D.dep+P.hro)).^2,1)); % x coord.(1st) [m]
y2q = linspace(yc+D.dep+P.hri,yc,np);% y coord.(2nd) [m]
x2q = (P.rci-(D.dep+P.dwi)*sqrt(1-min(((y2q-yc)./(D.dep+P.hri)).^2,1))); % x coord.(2nd) [m]
y3q = linspace(-yc,-yc-D.dep-P.hri,np);% y coord.(3rd) [m]
x3q = flip(x2q);
y4q = linspace(-yc-D.dep-P.hro,-yc,np);% y coord.(4th) [m]
x4q = flip(x1q);
x = [P.ro,P.ro,...
    x1q,P.rco,P.rci,...
    x2q,P.ra,P.ra,...
    x3q,P.rci,P.rco,x4q];
y = [-P.lc/2,P.lc/2,y1q,...
    P.lc/2+P.hro+D.dep,P.lc/2+P.hri+D.dep,...
    y2q,P.lc/2,-P.lc/2,y3q,...
    -P.lc/2-P.hri-D.dep,-P.lc/2-P.hro-D.dep,y4q];
% plot conductors cross sectional view
h1 = fill3(x,zeros(size(x)),y,colormap([0.9531,0.3008,0.1719])); % plot conductors
hold on;
fill3(-x,zeros(size(x)),y,colormap([0.9531,0.3008,0.1719])); % plot conductors
% plot conductor 3D
q = linspace(0,pi,np);
[rr,qq] = meshgrid(x,q);
surf(rr.*cos(qq),rr.*sin(qq),repmat(y,[size(rr,1),1]),...
    'FaceColor',colormap([0.9531,0.3008,0.1719]),'EdgeColor','None');
% epoxy protective layer
q1 = linspace(0,pi/2,np)';      % theta 1st quadrant [rad]
q2 = linspace(pi/2,pi,np)';     % theta 2nd quadrant [rad]
q3 = linspace(pi,3*pi/2,np)';   % theta 3rd quadrant [rad]
q4 = linspace(3*pi/2,2*pi,np)'; % theta 4th quadrant [rad]
x1q = D.dep*cos(q1)+P.rco;      % x coordinates (1st quadrant) [m]
y1q = D.dep*sin(q1)+P.lc/2;     % y coordinates (1st quadrant) [m]
x2q = D.dep*cos(q2)+P.rci;      % x coordinates (2nd quadrant) [m]
y2q = D.dep*sin(q2)+P.lc/2;     % y coordinates (2nd quadrant) [m]
x3q = D.dep*cos(q3)+P.rci;      % x coordinates (3rd quadrant) [m]
y3q = D.dep*sin(q3)-P.lc/2;     % y coordinates (3rd quadrant) [m]
x4q = D.dep*cos(q4)+P.rco;      % x coordinates (4th quadrant) [m]
y4q = D.dep*sin(q4)-P.lc/2;     % y coordinates (4th quadrant) [m]
x = [P.rwo;P.rwo;...
    x1q;P.rco;P.rci;...
    x2q;P.rei;P.rei;...
    x3q;P.rci;P.rco;x4q];
y = [-P.lc/2;P.lc/2;y1q;...
    P.lc/2+D.dep;P.lc/2+D.dep;...
    y2q;P.lc/2;-P.lc/2;y3q;...
    -P.lc/2-D.dep;-P.lc/2-D.dep;y4q];
h2 = fill3(x,zeros(size(x)),y,colormap([3, 86, 31]./255)); % plot protective layer
hold on;
fill3(-x,zeros(size(x)),y,colormap([3, 86, 31]./255)); % plot protective layer

% if Bpk is parsed plot relative mangetic permeability
if nargin>3
    cmap = [get(gcf,'Colormap');parula(256)];
    colormap(cmap);
    r = linspace(P.rci,P.rco,numel(Bpk));   % core radii for plotting over x
    l = linspace(-P.lc/2,P.lc/2,numel(Bpk));% core height for plotting over z
    [xx,zz] = meshgrid(r,l);        % create coordinate grid for plotting
    perf    = ppval(P.mur,r);       % fitted permeability function
    mur = transpose(repmat(fmuB(P.mp,Bpk,perf(:))/(4e-7*pi),[1,numel(Bpk)]));
    hold on;
    offset = P.ro/1e3;
    h3 = surf(xx,zeros(size(xx))-offset,zz,mur,'FaceColor','interp','EdgeColor','None');
    surf(-xx,zeros(size(xx))-offset,zz,mur,'FaceColor','interp','EdgeColor','None');
    hold off;
    cb = colorbar('location','eastoutside');
    caxis([min(min(mur));max(max(mur))]);
    cblim = get(cb,'Limits');
    set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
    title(cb,'$\mu_r$','interpreter','latex');
    cb.Ruler.TickLabelFormat = '%.2f';
else
    % plot core
    x = [ P.rci;  P.rco;   P.rco; P.rci];
    y = [-1;-1;1;1]*P.lc/2;
    h3 = fill3(x,zeros(size(x)),y,colormap([0.4336,0.5117,0.6289])); % plot core
    hold on;
    fill3(-x,zeros(size(x)),y,colormap([0.4336,0.5117,0.6289])); % plot core
    hold off;
end
grid on; grid minor;
title('3D cross section view, [m]','interpreter','latex');
legend([h1 h2 h3],{'Conductors','Plastic Protection','Core'},'interpreter','latex',...
    'location','best');
material metal;
lighting gouraud;
daspect([1 1 1]);
view(45,25);
xlabel('$x$, [m]','interpreter','latex');
ylabel('$y$, [m]','interpreter','latex');
zlabel('$z$, [m]','interpreter','latex');
camlight(0,0,'infinite');
camlight(0,0,'infinite');
camlight(0,0,'infinite');
end

%%
function fn = PlotTemperature(fn,P,D,T)

% define local variables
np = 100;

% plot cross section
set(figure(fn),'color','white'); fn=fn+1;
yce  = P.lc/2;
y1q = linspace(yce,yce+D.dep+P.hro,np);% y coord.(1st) [m]
x1q = P.rco+(D.dep+P.dwo)*sqrt(1-min(((y1q-yce)./(D.dep+P.hro)).^2,1)); % x coord.(1st) [m]
y2q = linspace(yce+D.dep+P.hri,yce,np);% y coord.(2nd) [m]
x2q = (P.rci-(D.dep+P.dwi)*sqrt(1-min(((y2q-yce)./(D.dep+P.hri)).^2,1))); % x coord.(2nd) [m]
y3q = linspace(-yce,-yce-D.dep-P.hri,np);% y coord.(3rd) [m]
x3q = flip(x2q);
y4q = linspace(-yce-D.dep-P.hro,-yce,np);% y coord.(4th) [m]
x4q = flip(x1q);
xcd = [P.ro,P.ro,...
    x1q,P.rco,P.rci,...
    x2q,P.ra,P.ra,...
    x3q,P.rci,P.rco,x4q];
ycd = [-P.lc/2,P.lc/2,y1q,...
    P.lc/2+P.hro+D.dep,P.lc/2+P.hri+D.dep,...
    y2q,P.lc/2,-P.lc/2,y3q,...
    -P.lc/2-P.hri-D.dep,-P.lc/2-P.hro-D.dep,y4q];
hold on;
% epoxy protective layer
q1 = linspace(0,pi/2,np)';      % theta 1st quadrant [rad]
q2 = linspace(pi/2,pi,np)';     % theta 2nd quadrant [rad]
q3 = linspace(pi,3*pi/2,np)';   % theta 3rd quadrant [rad]
q4 = linspace(3*pi/2,2*pi,np)'; % theta 4th quadrant [rad]
x1q = D.dep*cos(q1)+P.rco;      % x coordinates (1st quadrant) [m]
y1q = D.dep*sin(q1)+P.lc/2;     % y coordinates (1st quadrant) [m]
x2q = D.dep*cos(q2)+P.rci;      % x coordinates (2nd quadrant) [m]
y2q = D.dep*sin(q2)+P.lc/2;     % y coordinates (2nd quadrant) [m]
x3q = D.dep*cos(q3)+P.rci;      % x coordinates (3rd quadrant) [m]
y3q = D.dep*sin(q3)-P.lc/2;     % y coordinates (3rd quadrant) [m]
x4q = D.dep*cos(q4)+P.rco;      % x coordinates (4th quadrant) [m]
y4q = D.dep*sin(q4)-P.lc/2;     % y coordinates (4th quadrant) [m]
xep = [P.rwo;P.rwo;...
    x1q;P.rco;P.rci;...
    x2q;P.rei;P.rei;...
    x3q;P.rci;P.rco;x4q];
yep = [-P.lc/2;P.lc/2;y1q;...
    P.lc/2+D.dep;P.lc/2+D.dep;...
    y2q;P.lc/2;-P.lc/2;y3q;...
    -P.lc/2-D.dep;-P.lc/2-D.dep;y4q];
% plot region edge lines
xcr = [ P.rci;  P.rco;   P.rco; P.rci; P.rci];
ycr = [-1;-1;1;1;-1]*P.lc/2;
plot(xcd,ycd,'-k',xep,yep,'--k',xcr,ycr,'-.k','LineWidth',2); % plot core

%%
% plot core
for k = 1:D.ns
    % plot core region
    RegionTemperature(linspace(P.rsec(k),P.rsec(k+1),np),...
        linspace(0,P.lc/2,np),linspace(0,P.lc/2,np),T.core(k).Tcf);
    RegionTemperatureMirror(linspace(P.rsec(k),P.rsec(k+1),np),...
        linspace(0,P.lc/2,np),linspace(-P.lc/2,0,np),T.core(k).Tcf);
    % plot winding top region
    RegionTemperature(linspace(P.rsec(k),P.rsec(k+1),np),...
        linspace(0,P.hwm(k),np),linspace(D.dep+P.lc/2,...
        D.dep+P.lc/2+P.hwm(k),np),T.winding.top(k).Tcf);
    % plot winding bottom region
    RegionTemperature(linspace(P.rsec(k),P.rsec(k+1),np),...
        linspace(0,P.hwm(k),np),linspace(-D.dep-P.lc/2,...
        -D.dep-P.lc/2-P.hwm(k),np),T.winding.top(k).Tcf);
end
%%
% plot inner winding
RegionTemperature(linspace(P.ra,P.rei,np),linspace(0,P.lc/2,np),...
    linspace(0,P.lc/2,np),T.winding.ivert.Tcf);
RegionTemperatureMirror(linspace(P.ra,P.rei,np),linspace(0,P.lc/2,np),...
    linspace(-P.lc/2,0,np),T.winding.ivert.Tcf);
% plot outer winding
RegionTemperature(linspace(P.rwo,P.ro,np),linspace(0,P.lc/2,np),...
    linspace(0,P.lc/2,np),T.winding.overt.Tcf);
RegionTemperatureMirror(linspace(P.rwo,P.ro,np),linspace(0,P.lc/2,np),...
    linspace(-P.lc/2,0,np),T.winding.overt.Tcf);
% plot inner upper corner
le = (pi*P.dwi*P.rci/2-2*P.dwi^2/3)/(pi*P.rci-2*P.dwi);
rm = le*(pi*P.rci-2*P.dwi)/P.dwi;
dr = P.dwi^2/(4*le);
RegionTemperatureCorner(linspace(rm-dr,rm+dr,np),...
    linspace(0,le,np),linspace(D.dep,D.dep+P.dwi,np),...
    linspace(pi/2,pi,np),P.rci,P.lc/2,T.winding.ic.Tcf);
% plot inner lower corner
RegionTemperatureCorner(linspace(rm-dr,rm+dr,np),...
    linspace(0,le,np),linspace(D.dep,D.dep+P.dwi,np),...
    linspace(pi,3*pi/2,np),P.rci,-P.lc/2,T.winding.ic.Tcf);
% plot outer upper corner
le = (pi*P.dwo*P.rco/2+2*P.dwo^2/3)/(pi*P.rco+2*P.dwo);
rm = le*(pi*P.rco+2*P.dwo)/P.dwo;
dr = P.dwo^2/(4*le);
RegionTemperatureCorner(linspace(rm-dr,rm+dr,np),...
    linspace(0,le,np),linspace(D.dep,D.dep+P.dwo,np),...
    linspace(0,pi/2,np),P.rco,P.lc/2,T.winding.oc.Tcf);
% plot outer lower corner
RegionTemperatureCorner(linspace(rm-dr,rm+dr,np),...
    linspace(0,le,np),linspace(D.dep,D.dep+P.dwo,np),...
    linspace(3*pi/2,2*pi,np),P.rco,-P.lc/2,T.winding.oc.Tcf);
% plot symmetry axis
plot(zeros(2,1),[-1,1]*(D.dep+P.lc/2+P.dwi)*1.2,'-.r');
hold off;
% camlight;
cb = colorbar('location','eastoutside');
cblim = get(cb,'Limits');
caxis([cblim(1),cblim(2)]);
set(cb,'Ticks',linspace(cblim(1),cblim(2),9));
cb.Ruler.TickLabelFormat = '%.2f';
colormap('hot');
title(cb,'$T$ [$^{\circ}C$]','interpreter','latex');
view(0,90);
shading interp;
lighting gouraud;
xlim([0,(P.rco+P.dwo)*1.2]);
ylim([-1,1]*(P.lc/2+P.dwi)*1.2);
daspect([1 1 1]);
xlabel('$x$, [m]','interpreter','latex');
ylabel('$y$, [m]','interpreter','latex');
zlabel('$z$, [m]','interpreter','latex');
title(['Axisymmetric cross section view - $T_{amb}=$ ',...
    num2str(D.Ta-273.15),' [$^{\circ}C$]'],'interpreter','latex');
legend({'Winding','Epoxy Layer','Core'},'interpreter','latex',...
    'location','best');

end

%%
function [] = RegionTemperature(r,zcalc,zplot,Tcf)
% compute temperature profile in the region
[rr,zc] = meshgrid(r,zcalc);
Treg = (Tcf.c2r*rr.^2+Tcf.clr*log(rr)+Tcf.c2z*zc.^2+Tcf.c1z*zc+...
    Tcf.c0);
% plot temperature profile in the region
[rr,zp] = meshgrid(r,zplot);
surf(rr,zp,Treg-273.15,'EdgeColor','None');
end

%%
function [] = RegionTemperatureCorner(rcalc,zcalc,rplot,qplot,xc,yc,Tcf)
% compute temperature profile in the region
[rc,zc] = meshgrid(rcalc,zcalc);
Treg = transpose(Tcf.c2r*rc.^2+Tcf.clr*log(rc)+Tcf.c2z*zc.^2+Tcf.c1z*zc+...
    Tcf.c0);
% plot temperature profile in the region
[rp,qp] = meshgrid(rplot,qplot);
surf(rp.*cos(qp)+xc,rp.*sin(qp)+yc,Treg-273.15,'EdgeColor','None');
end

%%
function [] = RegionTemperatureMirror(r,zcalc,zplot,Tcf)
% compute temperature profile in the region
[rr,zc] = meshgrid(r,zcalc);
Treg = (Tcf.c2r*rr.^2+Tcf.clr*log(rr)+Tcf.c2z*zc.^2+Tcf.c1z*zc+...
    Tcf.c0);
% plot temperature profile in the region
[rr,zp] = meshgrid(r,zplot);
surf(rr,zp,flipud(Treg)-273.15,'EdgeColor','None');
end % end of RegionTemperatureMirror

%%
function [] = CheckDesignSpecsD(D)

if ~isfield(D,'Lrqi')
    error('Please, define D.Lrqi: minimum required incremental inductance [H].');
elseif ~isfield(D,'kpfmx')
    error('Please, define D.kpfmx: maximum allowed packing factor.');
elseif ~isfield(D,'kb')
    error('Please, define D.kb: winding build factor, determines distance between two consecutive conductors.');
elseif ~isfield(D,'amxa')
    error('Please, define D.amxa: maximum aspect ratio.');
elseif ~isfield(D,'MLmxa')
    error('Please, define D.MLmxa: maximum allowed mass [kg].');
elseif ~isfield(D,'PLmxa')
    error('Please, define D.MLmxa: maximum allowed power dissipation [W].');   
elseif ~isfield(D,'Bmxa')
    error('Please, define D.Bmxa: maximum allowed flux density [T].');
elseif ~isfield(D,'nspc')
    error('Please, define D.nspc: # of parallel strands per conductor.');  
elseif ~isfield(D,'dep')
    error('Please, define D.dep: core protection layer thickness [m].');  
elseif ~isfield(D,'lwrmx')
    error('Please, define D.lwrmx: maximum allowed wire length [m].');  
elseif ~isfield(D,'dsmx')
    error('Please, define D.dsmx: maximum allowed strand diameter [m].');  
elseif ~isfield(D,'ep')
    if ~isfield(D.ep,'rho')
        error('Please, define D.ep.rho: protection layer mass density [kg/m^3].'); 
    end
end

end % end of CheckDesignSpecs()