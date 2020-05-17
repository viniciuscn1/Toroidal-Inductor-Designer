function [T,Pw,Pp,c,Rdc] = ToroidTEC(D,P,Pls)
%%
% define local constants
mu0     = 4e-7*pi;              % permeability of free space [H/m]
P.w     = D.wn;                 % harmonics radian frequencies [rad/s]
P.irmsn = D.irmsn;              % harmonics rms currents [A]

% initialize local variables
c = 1;                          % convergence flag

%%
% initialize TEC
Nsecnodes = 5*D.ns+1;           % # of sectional nodes (depends on # of
                                % sections): 5*D.ns+1 core
Nsecwtop  = 6*D.ns+1;           % 6*D.ns+1 winding top sections
Nnodes_wv = 6;                  % # of nodes vertical winding regions
Nnodes_corners = 4;             % # of nodes winding corners
Nnodes = Nsecnodes+Nsecwtop+(Nnodes_wv+Nnodes_corners)*2;
TEC = tec_init(Nnodes);         % initialize TEC structure

% add material for the core
[TEC,micr] = tec_material(TEC,'Core',P.mp.TEC.kx,P.mp.TEC.ky,...
    P.mp.TEC.kz,TEC.ml.c(3),P.mp.rho);

% look up material properties
ka = TEC.ml.kx(D.mia);          % thermal conductivity of air [W/(m.K)]
% compute heat transfer coefficient of core to winding
hcw = 4*D.hsl*ka/(4*ka+D.hsl*(4*D.gwc+(4-pi)*P.rs));

% compute homogenized winding thermal conductivities ----------------------
% inner side vertical section of the winding
[TEC,miwiv] = tec_homogenize_rwb(TEC,pi*(P.ra+P.rei),P.dwi,P.N*D.nspc,...
                P.rs,D.ti,D.micd,D.mii,D.mia,'inner vertical');
                         
% inner corner thermal properties
miwic            = TEC.ml.nm+1;     % increment material index
TEC.ml.nm        = miwic;           % assign material index
TEC.ml.kx(miwic) = TEC.ml.kz(miwiv);% inner vertical kz->inner corner kx
TEC.ml.ky(miwic) = TEC.ml.kz(miwiv);% inner vertical kz->inner corner ky
TEC.ml.kz(miwic) = TEC.ml.kx(miwiv);% inner vertical kx->inner corner kz
TEC.ml.c(miwic)  = TEC.ml.c(miwiv); % inner vertical c->inner corner c
TEC.ml.row(miwic)= TEC.ml.row(miwiv);% inner vertical mass density ->inner corner mass density
TEC.ml.des{miwic}= 'inner corner';

% outer side vertical section of the winding
[TEC,miwov] = tec_homogenize_rwb(TEC,pi*(P.rwo+P.ro),P.dwo,P.N*D.nspc,...
                P.rs,D.ti,D.micd,D.mii,D.mia,'outer vertical');
                         
% outer corner thermal properties
miwoc            = TEC.ml.nm+1;     % increment material index
TEC.ml.nm        = miwoc;           % assign outer corner material index
TEC.ml.kx(miwoc) = TEC.ml.kz(miwov);% outer vertical kz->outer corner kx
TEC.ml.ky(miwoc) = TEC.ml.kz(miwov);% outer vertical kz->outer corner ky
TEC.ml.kz(miwoc) = TEC.ml.kx(miwov);% outer vertical kx->outer corner kz
TEC.ml.c(miwoc)  = TEC.ml.c(miwov); % outer vertical c-> outer corner c
TEC.ml.row(miwoc)= TEC.ml.row(miwov);% outer vertical mass density->outer corner mass density
TEC.ml.des{miwoc}= 'outer corner';

% top and bottom sections of the winding
miwtb = zeros(D.ns,1);  % initialize material idx winding top and bottom
for k = 1:D.ns          % loop through however many pre-defined sections
    [TEC,miwtb(k)] = tec_homogenize_rwb(TEC,2*pi*P.rsecm(k),P.hwm(k),...
        P.N*D.nspc,P.rs,D.ti,D.micd,D.mii,D.mia,'homogenized winding');
    % since the conductors are oriented in the radial direction in this
    % section, it is necessary to flip the thermal conductivity of radial
    % and axial orientations
    TEC.ml.kx(miwtb(k)) = TEC.ml.kz(miwtb(k));
    TEC.ml.kz(miwtb(k)) = TEC.ml.ky(miwtb(k));
    TEC.ml.ky(miwtb(k)) = TEC.ml.kx(miwtb(k));
end

%%
% setup TEC

% core region -------------------------------------------------------------
C   = cell(D.ns,1);     % initialize each core section parameters struct
Wt  = cell(D.ns,1);     % initialize each winding top sections struct
for k = 1:D.ns          % loop through each section
    % define each cylindrical section of the core as a TEC region
    [TEC,C{k}] = tec_cylindrical(TEC,micr,P.rsec(k),P.rsec(k+1),P.lc/2,2*pi,...
        5*k-4,-1,5*k-2,5*k-1,5*k,5*k+1,5*k-3,Pls(k)/2,0);
    % define each cylindrical section of the winding top as a TEC region
    [TEC,Wt{k}] = tec_cylindrical(TEC,miwtb(k),P.rsec(k),P.rsec(k+1),P.hwm(k),...
        2*pi,Nsecnodes+6*k-5,Nsecnodes+6*k-4,Nsecnodes+6*k-2,...
        Nsecnodes+6*k-1,Nsecnodes+6*k,Nsecnodes+6*k+1,Nsecnodes+6*k-3,0,0);
    % define upper region of each winding section in the toroid top
    TEC = tec_ambient_branch(TEC,Nsecnodes+6*k-3,Wt{k}.Sz,D.hwa,D.Ta);
    % define contact region between core and top winding section
    TEC = tec_contact_branch(TEC,Nsecnodes+6*k-4,5*k-3,C{k}.Sz,hcw);
end

% inner vertical winding region -------------------------------------------
nividx = Nsecnodes+Nsecwtop+1;  % starting node inner vertical winding reg.
[TEC,Wiv] = tec_cylindrical(TEC,miwiv,P.ra,P.rei,P.lc/2,2*pi,...
        nividx,-1,nividx+2,nividx+3,nividx+4,nividx+5,nividx+1,...
        0,0);
% inner vertical winding region to core contact thermal resistance
TEC = tec_contact_branch(TEC,1,nividx+5,Wiv.Sro,hcw);
% inner vertical winding region to ambient
TEC = tec_ambient_branch(TEC,nividx,Wiv.Sri,D.hwa,D.Ta);

% outer vertical winding region -------------------------------------------
novidx = nividx+Nnodes_wv;      % starting node outer vertical winding reg.
[TEC,Wov] = tec_cylindrical(TEC,miwov,P.rwo,P.ro,P.lc/2,2*pi,...
        novidx,-1,novidx+2,novidx+3,novidx+4,novidx+5,novidx+1,...
        0,0);
% outer vertical winding region to core contact thermal resistance
TEC = tec_contact_branch(TEC,5*D.ns+1,novidx,Wov.Sri,hcw);
% outer vertical winding region to ambient
TEC = tec_ambient_branch(TEC,novidx+5,Wov.Sro,D.hwa,D.Ta);

% inner corner winding region ---------------------------------------
nicidx = novidx+Nnodes_wv;
wic = (P.dwi+P.hri)/2;
le = (pi*wic*P.rci/2-2*wic^2/3)/(pi*P.rci-2*wic);
rm = le*(pi*P.rci-2*wic)/wic;
dr = wic^2/(4*le);
[TEC,Wic] = tec_cylindrical(TEC,miwic,rm-dr,rm+dr,le,2*pi,...
        nividx+1,-1,nicidx+1,nicidx+2,nicidx+3,...
        Nsecnodes+1,nicidx,0,0);
TEC = tec_ambient_branch(TEC,nicidx,Wic.Sz,D.hwa,D.Ta);  % ir to Tamb

% outer corner winding region ---------------------------------------
nocidx = nicidx+Nnodes_corners;
woc = (P.dwo+P.hro)/2;
le = (pi*woc*P.rco/2+2*woc^2/3)/(pi*P.rco+2*woc);
rm = le*(pi*P.rco+2*woc)/woc;
dr = woc^2/(4*le);
[TEC,Woc] = tec_cylindrical(TEC,miwoc,rm-dr,rm+dr,le,2*pi,...
        novidx+1,-1,nocidx+1,nocidx+2,nocidx+3,...
        Nsecnodes+6*D.ns+1,nocidx,0,0);
TEC = tec_ambient_branch(TEC,nocidx,Woc.Sz,D.hwa,D.Ta);  % ir to Tamb

%%
% solve TEC

% In this problem the rms current is fixed, hence Jrms = Irms/ac. As the
% coil is heated, its resistance increases. For a fixed current, the
% resistive power loss also increases. This requires some iterations of the
% solver to reach steady state.

% store initial TEC
TEC0 = TEC;

% initialize temperature in all nodes as the ambient temperature
Tsol = ones(TEC.nn,1)*D.Ta;     % initialize solution nodes as ambient temp.
Told = Tsol;                    % save temperature solution array 
T0 = P.cp.t0+273.15;            % ambient temperature [K]
didt2 = P.w.^2.*D.irmsn.^2;     % time average of (didt)^2

%%
% iteratively solve for temperature
k  = 1;                 % initialize iterations counter
re = 0;                 % initialize residual error as zero
Rhts = zeros(D.ns,1);   % initialize ohmic losses on top winding sections
Phts = zeros(D.ns,1);   % initialize ohmic losses on top winding sections
Pphts = zeros(D.ns,1);  % initialize proximity losses top winding sections
while (k==1)||((k<=D.kmax)&&(re>D.emax))
    TEC = TEC0;         % restore initial TEC
    
    % compute updated conductivity for each region
    sigmaiv  = P.cp.sigma0/(1+P.cp.alpha*(Tsol(nividx+4)-T0));
    sigmaov  = P.cp.sigma0/(1+P.cp.alpha*(Tsol(novidx+4)-T0));
    sigmaic  = P.cp.sigma0/(1+P.cp.alpha*(Tsol(nicidx+3)-T0));
    sigmaoc  = P.cp.sigma0/(1+P.cp.alpha*(Tsol(nocidx+3)-T0));
    
    % compute resistive losses in every section of the coil (ohmic+skin
    % effect)
    [Piv,Riv] = WindingLosses(P,sigmaiv,P.Vcdi,D.nspc);    % inner vertical
    [Pov,Rov] = WindingLosses(P,sigmaov,P.Vcdo,D.nspc);    % outer vertical
    [Pic,Ric] = WindingLosses(P,sigmaic,P.Vcdic,D.nspc);   % inner corner
    [Poc,Roc] = WindingLosses(P,sigmaoc,P.Vcdoc,D.nspc);   % outer corner

    % compute proximity losses
    % compute the dynamic resistances
    Dov = mu0*P.N^3*pi*sigmaov*P.rs^4*P.lc*P.Pe_ov_int/(4*P.Vclo);
    Div = mu0*P.N^3*pi*sigmaiv*P.rs^4*P.lc*P.Pe_iv_int/(4*P.Vcli);
    
    % compute proximity losses [W]
    Ppiv = sum(Div*didt2);	% sum all harmonic components of proximity losses 
    Ppov = sum(Dov*didt2);  % sum all harmonic components of proximity losses 

    % update power dissipation in every section of the coil in the TEC
    TEC.Pn(nividx+4) = (Piv+Ppiv)/2; % inner vertical region (half for symmetry)
    TEC.Pn(novidx+4) = (Pov+Ppov)/2; % outer vertical region (half for symmetry)
    TEC.Pn(nicidx+3) = Pic;          % inner corner
    TEC.Pn(nocidx+3) = Poc;          % outer corner
    
    % update sectionalized regions of the winding
    for i = 1:D.ns          % loop through each section         
        % update conductivity of horizontal top section
        sigmaht = P.cp.sigma0/(1+P.cp.alpha*(Tsol(Nsecnodes+6*i)-T0));
        % compute proximity losses for top winding section
        Dht = mu0*P.N^3*pi*sigmaht*P.rs^4*P.dc/...
            D.ns.*P.Pe_hs(i)./(4*P.Vcls(i));
        Pphts(i) = sum(Dht*didt2);
        % compute ohmic and skin effect losses for top winding section
        [Phts(i),Rhts(i)] = WindingLosses(P,sigmaht,P.Vcds(i),D.nspc);
        % update power dissipated in top winding section
        TEC.Pn(Nsecnodes+6*i) = Phts(i)+Pphts(i);    % top sections
    end
    
    % solve TEC
    Tsol = tec_linear_solve(TEC);   % call TEC linear solver
    
    % identify unstable numerical problem
    if      any((Tsol-D.Ta)<0)||...
            ((min(Tsol)<0)&&(max(abs(Told-Tsol))>re))||...
            ((re>1e5)&&(max(abs(Told-Tsol))>re))
        c       = 0;
        Rdc     = 0;
        Pw      = inf;
        Pp      = inf;
        T.Twmx  = inf;
        T.Tcrmx = inf;
        return;
    end
    re   = max(abs(Told-Tsol));     % compute residual error
    Told = Tsol;                    % save current solution as old
    k    = k+1;                     % increment counter
end

%%
% post-processing 

% compute final value of proximity effect losses
Pp = Ppiv+Ppov+sum(Pphts)*2;

% initialize mean and peak temperature output variables
T.core(D.ns).Tmn = 0;
T.winding.top(D.ns).Tmn = 0;
Tcrpk = zeros(D.ns,1);                      % core peak temperature [K]
Twpk = zeros(4+D.ns,1);                     % winding peak temperature [K]
% get mean and peak temperature in the core and winding
for k = 1:D.ns                              % loop through each section
    [T.core(k).Tmn,T.core(k).Tpk,...
        T.core(k).Tcf] = ...                % get mean and peak temperature    
        tec_cylindrical_temp(C{k},Tsol);    % of core sections
    [T.winding.top(k).Tmn,T.winding.top(k).Tpk,...
        T.winding.top(k).Tcf] = ...  
        tec_cylindrical_temp(Wt{k},Tsol);   % of winding top sections
    Tcrpk(k) = T.core(k).Tpk;               % store core section peak temp.
    Twpk(k)  = T.winding.top(k).Tpk;        % store winding section peak temp.
end
[T.winding.ivert.Tmn,T.winding.ivert.Tpk,...% winding inner vertical section
    T.winding.ivert.Tcf] = tec_cylindrical_temp(Wiv,Tsol);    
[T.winding.overt.Tmn,T.winding.overt.Tpk,...% winding outer vertical section
    T.winding.overt.Tcf] = tec_cylindrical_temp(Wov,Tsol);     
[T.winding.ic.Tmn,T.winding.ic.Tpk,...      % winding inner corner section
    T.winding.ic.Tcf] = tec_cylindrical_temp(Wic,Tsol);     
[T.winding.oc.Tmn,T.winding.oc.Tpk,...    % winding outer corner section
    T.winding.oc.Tcf] = tec_cylindrical_temp(Woc,Tsol);    
% store winding corner and vertical sections peak temperatures
Twpk(end-3:end) = [T.winding.ivert.Tpk;T.winding.overt.Tpk;T.winding.ic.Tpk;
    T.winding.oc.Tpk];

% find maximum temperature
Twmx = max(Twpk);               % winding maximum temperature [K]
Tcrmx= max(Tcrpk);              % core maximum temperature [K]

% determine winding ohmic and skin effect losses [W]
Pw  = Piv+Pov+(Pic+Poc+sum(Phts))*2;
Rdc = Riv+Rov+(Ric+Roc+sum(Rhts))*2;

%%
% assign output variables to output structure
T.Twpk  = Twpk;             % peak temperatures in winding regions [K]
T.Tcrpk = Tcrpk;            % peak temperatures in core sections [K]
T.Tsol  = Tsol;             % temperature in every node [K]
T.Tmx   = max(Tcrmx,Twmx);  % maximum device temperature [K]
T.Twmx  = Twmx;             % maximum winding temperature [K]
T.Tcrmx = Tcrmx;            % maximum core temperature [K]

% check for convergence
if isnan(Twmx)
    c = 0;
    Pw = inf;
    Pp = inf;
    Rdc= 0;
    T.Twmx = inf;
    T.Tcrmx = inf;
    return;
end

end

%%
% compute winding losses (ohmic+skin effect)
function [Pw,Rdc] = WindingLosses(P,sigma,Vcd,nspc)
%% WindingLosses.m
% This function computes the winding ohmic losses. This losses include the
% losses due to the traditional fundamental and its harmonic components,
% thus accounting for additional losses from skin effect. The greater the
% harmonic frequency, the larger the apparent AC resistance of the wire is
% due to the skin depth, which results in greater losses than if the entire
% cross section of wire had been considered.
%
% Call: 
%  Pac = WindingLosses(P,W,sigma,Vcd);
% 
% Input:
%  P        = toroidal core parameter structure
%   P.w     = harmonic radian frequencies [rad/s]
%   P.irmsn = rms current of each harmonic component [A]
%  W        = winding parameter structure
%   P.rc    = conductor radius [m]
%   P.ac    = conductor cross sectional area [m^2]
%  sigma    = wire conductivity [S/m]
%  Vcd      = conductor volume [m^3]
% 
% Output:
%  Pac      =  AC ohmic power losses [W]
% 
% Internal Variables:
%  mu0      = permeability of free space [H/m]
%  kappa    = defined variable for Bessel function
%  rhat     = transformed radius for Bessel function
%  drhat    = used for numeric differentiation (1st order finite difference)
%  J0       = Bessel function solution for rhat with order v=0
%  J02      = Bessel function solution for rhat+drhat with order v=0
%  J0p      = finite difference of Bessel function 
%  Z        = AC impedance [Ohm]
%  Rac      = AC resistance [Ohm]
%  Pacn     = AC ohmic losses for each individual harmonic component [W]
% 
% Written by: Vinicius Nascimento
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(P.w~=0)  % AC    
    % compute skin effect losses parameters
    kappa = sqrt(1i./(P.w(P.w~=0)*sigma*4e-7*pi)); % defined variable for Bessel function
    rhat = P.rs./kappa;         % transformed radius for Bessel function
    drhat = rhat*1e-6;          % used for numeric differention
    J0 = besselj(0,rhat);       % solve Bessel function for rhat with order v=0
    J02= besselj(0,rhat+drhat); % solve Bessel function for slightly greater rhat
    J0p= (J02-J0)./drhat;       % solve finite difference
    Z = -(Vcd/P.ac)*J0./(2*pi*P.rs*sigma*kappa.*J0p); % compute ac impedance
    Rwac = real(Z);             % get strand AC resistance [Ohms]
    Rac = Rwac/nspc;            % total coil AC resistance [Ohms]
    
    % compute ac losses (including ohmic + skin effect) [W]
    Pacn = Rac.*(P.irmsn(P.w~=0)).^2;   % ac losses for each harmonic [W]
    Rdc = Vcd/(sigma*P.ac^2);
    Pdc = Rdc*P.irmsn(P.w==0)^2;
    Pac = sum(Pacn);
    % compute total ac losses (ohmic + skin effect) [W]
    if ~isempty(Pdc)
        Pw = Pac+Pdc;
    else
        Pw = Pac;
    end
else    % DC
    Rdc = Vcd/(sigma*P.ac^2);
    Pw = Rdc*P.irmsn(P.w==0)^2;
end

end

%%
% updated tec_cylindrical.m function 
function [TEC,CE] = tec_cylindrical(TEC,mc,ri,ro,lz,theta, ...
                                    nir,n0z,ncr,ncz,nmn, ...
                                    nor,nlz,p,gp)                                                
% tec_cylindrical  Adds a hollow cylindrical element 
%                  to a thermal equivalent circuit
%
% [TEC]      = tec_cylindrical(TEC,mc,ri,ro,lz,theta,nir,n0z,ncr,ncz, ...
%                              nmn,nor,nlz,p,gp)
%
% Inputs: 
% TEC      = TEC data structure (see tec_init for documentation)
% mc       = material code
%            Note: thermal conuctivity in x direction is used as radial
%                  direction
% ri       = inner radius of cylinder (m)
% ro       = outer radius of cylinder (m)
% lz       = length in z-direction (m)
% theta    = angle spanned by cyllinder (m)
% nir      = node at r=ri (use -1 for open circuit)
% n0z      = node number at z=0 (use -1 for open circuit)
% ncr      = node number at center of radial axis of 'T' circuit 
% ncz      = node number at center of z-axis 'T' circuit
% nmn      = node number for mean temperature of cuboid 
% nor      = node at r=ro (use -1 for open circuit)
% nlz      = node number at z=lz (use -1 for open circuit)
% p        = nominal power into cuboid (W)
% gp       = slope of power dissipation P.r.t. mean temperature (W/K)
%
% Outputs:
% CE       = structure of cyclindrical element parameters
%  CE.V    = volume (m^3)
%  CE.M    = mass (kg)
%  CE.C    = heat capacity (J/K)
%  CE.Gx   = thermal conductance from edge to center in x-direction (W/K)
%  CE.Gy   = thermal conductance from edge to center in y-direction (W/K)
%  CE.Gz   = thermal conductance from edge to center in z-direction (W/K)
%  CE.Sx   = area perpendicular to x-axis (m^2)
%  CE.Sy   = area perpendicular to y-axis (m^2)
%  CE.Sz   = area perpendicular to z-axis (m^2)
% TEC      = Updated TEC data structure
%
% Written by: S.D. Sudhoff
% Modified by: Vinicius Nascimento

% compute volume, mass, thermal capacitance,dimensions
CE.Sri=theta*ri*lz;
CE.Sro=theta*ro*lz;
CE.Sz=0.5*theta*(ro^2-ri^2);
CE.V=CE.Sz*lz;
CE.M=CE.V*TEC.ml.row(mc);
CE.C=CE.M*TEC.ml.c(mc);
CE.ri=ri;
CE.ro=ro;
CE.lz=lz;
CE.theta=theta;

% record node numbers
CE.nir=nir;
CE.ncr=ncr;
CE.nor=nor;
CE.n0z=n0z;
CE.ncz=ncz;
CE.nlz=nlz;
CE.nmn=nmn;

% check kr=kx=ky
if TEC.ml.kx(mc)~=TEC.ml.ky(mc)
    error('kr = kx = ky, but kx = %.3f and ky = %.3f',...
        TEC.ml.kx(mc),TEC.ml.ky(mc))
end

% record thermal conductivities
CE.kr = TEC.ml.kx(mc);
CE.kz = TEC.ml.kz(mc);

% compute thermal conductances
ro2     = ro^2;
ri2     = ri^2;
ro2mri2 = ro2-ri2;
ro2pri2 = ro2+ri2;
lnrori  = log(ro/ri);
ro2ri2  = ro2*ri2;
twokrthetal = 2.0*CE.kr*theta*lz;
CE.Gir  = twokrthetal/(2*ro2*lnrori/ro2mri2-1);
CE.Gor  = twokrthetal/(1-2*ri2*lnrori/ro2mri2);
CE.Gtr  = -2*twokrthetal*ro2mri2/ ...
        (ro2pri2-4*ro2ri2*lnrori/ro2mri2);
CE.Gz   = 2*CE.kz*CE.Sz/lz;

% put in branches
if (nir>=0)TEC=tec_g_branch(TEC,nir,ncr,CE.Gir); end;
if (nor>=0)TEC=tec_g_branch(TEC,nor,ncr,CE.Gor); end;
if (n0z>=0)TEC=tec_g_branch(TEC,n0z,ncz,CE.Gz); end;
if (nlz>=0)TEC=tec_g_branch(TEC,nlz,ncz,CE.Gz); end;
TEC=tec_g_branch(TEC,ncr,nmn,CE.Gtr);
TEC=tec_g_branch(TEC,ncz,nmn,(-3*CE.Gz));
TEC=tec_std_branch(TEC,nmn,0,0,gp,nmn,0,p);

end

%%
% update tec_cylindrical_temp.m function
function [Tmn,Tpk,Tcf]=tec_cylindrical_temp(C,T)

% tec_cylindrical_temp calculates the steady-state spatial mean and peak
% temperature within a cylindrical region
%
% [TAmn,TApk] = tec_cylindrical_temp(C,T)
%
% Inputs:
% C        = structure of cylindrical element parameters 
%            see tec_cylindrical.m for documentation
% T        = vector of nodal temperatures (K)
%
% Outputs:
% TAmn     = spatial mean temperature of cylindrical element (K)
% TApk     = spatial peak temperature of cylindrical element (K)
%
% Internal:
% Qir      = heat transfer rate into r=ri surface (W)
% Q0z      = heat transfer rate into z=0 surface (W)
% Qor      = heat transfer rate out of r=ro surface (W)
% Qlz      = heat transfer rate out of z=lz surface (W)
% c2r      = radial spatial temperature coefficient (r^2 term) (K/m^2)
% clr      = radial spatial temperature coefficeint (ln term) (K)
% c1z      = z-axis spatial temperature coefficient 1 (K/m)
% c2z      = z-axis spatial temperature coefficient 2(K/m^2)
% c0       = spatial temperature coefficient 0 (K)
% ro2      = intermediate variable (m^2)
% ri2      = intermediate variable (m^2)
% lz2      = intermediate variable (m^2)
% re       = radial extremum temperature point condidate (m)
% ze       = z-axis extremum temperature point condidate (m)
% Tri      = intermediate variable (K)
% Tro      = intermediate variable (K)
% Ter      = extremum temperature component in radial direction (K)
% Tez      = extremem temperature component in z-axis (K)
%
% Written by: S.D. Sudhoff                               
% Modified by: Vinicius Nascimento                         
  
% compute mean temperature
Tmn=T(C.nmn);
   
% compute heat flows    
if (C.nir>0)
   Qir=(T(C.nir)-T(C.ncr))*C.Gir;
else
   if (C.nir==0)    
      Qir=-T(C.ncr)*C.Gir;
   else
      Qir=0;
   end
end
if (C.n0z>0)
   Q0z=(T(C.n0z)-T(C.ncz))*C.Gz;
else
   if (C.n0z==0)
      Q0z=-T(C.ncz)*C.Gz;
   else
      Q0z=0;
   end
end
if (C.nor>0)
   Qor=(T(C.ncr)-T(C.nor))*C.Gor;
else
   if (C.nor==0)
      Qor=T(C.ncr)*C.Gor;
   else
      Qor=0;
   end
end
if (C.nlz>0)
   Qlz=(T(C.ncz)-T(C.nlz))*C.Gz;
else
   if (C.nlz==0) 
      Qlz=T(C.ncz)*C.Gz;
   else
      Qlz=0;
   end
end

% compute spatial distribution terms
ro2=C.ro^2;
ri2=C.ri^2;
lz2=C.lz^2;
clr=(Qor*ri2-Qir*ro2)/(2*C.kr*C.V);
c1z=-Q0z./(C.kz*C.Sz);
c2r=(Qir-Qor)/(4*C.kr*C.V);
c2z=(Q0z-Qlz)/(2*C.kz*C.V);
c0=Tmn-(0.5*c2r*(ro2+ri2)+ ...
        clr*((ro2*log(C.ro)-ri2*log(C.ri))/(ro2-ri2)-0.5)+...
        c2z*lz2/3+c1z*lz2/2);
   
% compute extremum temperature for radial term
Tri=c2r*ri2+clr*log(C.ri);
Tro=c2r*ro2+clr*log(C.ro);
Ter=max(Tri,Tro);
if (clr*c2r<0)
   re=sqrt(-0.5*clr/c2r);
   if (re>C.ri)&&(re<C.ro)
      Ter=max([Tri Tro c2r*re^2+clr*log(re)]);
   end
end
   
% compute extremum temperature for z-axis
Talz=c2z*lz2+c1z*C.lz;
Tez=max(0,Talz);
if (c2z~=0)
   ze=-c1z/(2*c2z);
   if (ze>0)&&(ze<C.lz)
      Tez=max([0 Talz -0.25*c1z^2/c2z]);
   end
end
   
% compute the peak temperature in the sample
Tpk=Ter+Tez+c0;

% return temperature profile coefficients
Tcf.c0  = c0;
Tcf.clr = clr;
Tcf.c2r = c2r;
Tcf.c1z = c1z;
Tcf.c2z = c2z;

end
