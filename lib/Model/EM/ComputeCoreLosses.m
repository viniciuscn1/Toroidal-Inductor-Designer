function P = ComputeCoreLosses(P,D)
% compute core losses
if strcmp(D.cl.wave,'sine')         % sine wave approximation for core losses
    if D.cl.fcl~=0 % AC
        % define AC operating condition
        wcl = 2*pi*D.cl.fcl;        % radian frequency for core losses [rad/s]
        qf= 2*pi*D.cl.symmetry;     % final angle for analysis [rad]
        qe= linspace(0,qf,D.cl.ncl);% angles - electrical cycle [rad]
        t = qe/wcl;                 % time vector - electrical cycle [s]
        I = D.cl.Icl*sin(wcl*t);    % generate quarter cycle of current [A]
        B = SolveB(P,I,P.rsec);     % solve for flux density [T] (row->space, col->time)
        Bsec = (B(1:end-1,:)+B(2:end,:))/2;  % section B assumed to be average
        % core loss density in each section [W/m^3]
        if strcmp(D.cl.method,'MSE')
            pld = ComputeMSE(t,Bsec,D.cl.fcl,P.mp,D.cl.symmetry);% spatial average of B for each section [T]    
        elseif strcmp(D.cl.method,'iGSE')
            pld = ComputeiGSE(t,Bsec,D.cl.fcl,P.mp,D.cl.symmetry);
        else
            error('The core loss method %s is invalid.',D.cl.method);
        end
        P.Plcsec = P.Vsec.*pld;     % core losses for each section [W]
    else % DC
        P.Plcsec = zeros(size(P.Vsec));
    end
elseif strcmp(D.cl.wave,'full')     % full wave is processed
    if D.cl.fcl~=0 % AC
        t = linspace(D.t(1),D.t(end),D.cl.ncl);% create a time vector with ncl points
        t = t-t(1);                 % set initial time to zero
        I = interp1(D.t,D.Ic,t);    % interpolate ncl points out of current waveform [A]
        B = SolveB(P,I,P.rsec);     % solve for flux density [T] (row->space, col->time)
        Bsec = (B(1:end-1,:)+B(2:end,:))/2;  % section B assumed to be average
        % core loss density in each section [W/m^3]
        if strcmp(D.cl.method,'MSE')
            pld = ComputeMSE(t,Bsec,D.cl.fcl,P.mp,D.cl.symmetry);% spatial average of B for each section [T]
        elseif strcmp(D.cl.method,'iGSE')
            pld = ComputeiGSE(t,Bsec,D.cl.fcl,P.mp,D.cl.symmetry);
        else
            error('The core loss method %s is invalid.',D.cl.method);
        end
        P.Plcsec = P.Vsec.*pld;     % core losses for each section [W]
    else % DC
        P.Plcsec = zeros(size(P.Vsec));
    end
else
    error('Wave approximation ''%s'' is invalid.',D.cl.wave);
end

% total core losses [W]
P.Plc = sum(P.Plcsec);              

end % end of ComputeCoreLosses()

%%
function pld = ComputeiGSE(t,B,f,mp,sym)
% ComputeiGSE computes the core loss density of a material
%             Uses a variation of the improved generalized Steinmetz equation.

% define local constants
np = length(t);            % # of points in the waveform
% compute ki
ki = mp.SE.kh./((2.^(mp.SE.beta+1)).*pi.^(mp.SE.alpha-1).*(0.2761+1.7061./(mp.SE.alpha+1.354)));
% determine integral of the square of the time derivative of flux density
pB        = diff(B,1,2)./repmat(diff(t),[size(B,1),1]);
pB(:,np)  = pB(:,np-1);
int_dBdt2 = (1/sym)*trapz(t,pB.^2,2);
adBdt     = abs(pB);                    % absolute value of dBdt
% determine DeltaB
if sym<1 % symmetry is being used (e.g. 1/2-> positive half of a sine, 1/4-> positive quarter of a sine)
    DB = 2*max(B,[],2);
elseif sym==1 % symmetry is not being used, the full cycle is processed
    DB = max(B,[],2)-min(B,[],2); 
else
    error('The symmetry value (%.2f) should be between 0 and 1',sym);
end 
% compute power loss density
pld = mp.SE.ke*f*int_dBdt2; % initialize with eddy current losses
for k = 1:numel(mp.SE.kh)   % for each SE parameter term
    pld = pld+(1/sym)*f*ki(k)*trapz(t,(adBdt.^mp.SE.alpha(k)).*...
        repmat((DB.^(mp.SE.beta(k)-mp.SE.alpha(k))),[1,np]),2);
end
end % end of ComputeiGSE()

%%
function pld = ComputeMSE(t,B,f,mp,sym)
% ComputeMSE computes the core loss density of a material
%             Uses a variation of the modified Steinmetz equation.

% define local constants
np = length(t);             % # of points in the waveform

% determine integral of the square of the time derivative of flux density
pB        = diff(B,1,2)./repmat(diff(t),[size(B,1),1]);
pB(:,np)  = pB(:,np-1);
int_dBdt2 = (1/sym)*trapz(t,pB.^2,2);
% determine DeltaB
if sym<1 % symmetry is being used (e.g. 1/2-> positive half of a sine, 1/4-> positive quarter of a sine)
    DB = 2*max(B,[],2);
elseif sym==1 % symmetry is not being used, the full cycle is processed
    DB = max(B,[],2)-min(B,[],2); 
else
    error('The symmetry value (%.2f) should be between 0 and 1',sym);
end 
% determine equivalent frequency in MSE
feq = zeros(size(DB));
feq(DB~=0) = 2*int_dBdt2(DB~=0)./(DB(DB~=0)*pi).^2;
% compute power loss density
pld = mp.SE.ke*f*int_dBdt2; % initialize with eddy current losses
for k = 1:numel(mp.SE.kh)   % increment with each SE component
    pld = pld+mp.SE.kh(k)*(DB/2).^(mp.SE.beta(k)).*feq.^(mp.SE.alpha(k)-1)*f;
end
end % end of ComputeMSE()