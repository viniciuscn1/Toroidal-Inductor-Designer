function [varargout] = DefineCurrent(D,fn)
%% Coil excitation
% There are two options:
%  1) Define current waveform by setting a time vector (D.t) in [s] and 
%  a current vector (D.Ic) in [A].
%  2) Define an array of harmonic frequencies (D.fcf) in [Hz], harmonic 
%  currents (D.Icf) in [A] and another for phases (D.phf) in [rad], such 
%  that the current can be defined in time as follows:
%  ifit = sum(D.Icf.*cos(2*pi*D.fcf*D.t'+D.phf));

%%
% input pre-processing

% check if basic current description has been provided 
if (~isfield(D,'Ic'))||(~isfield(D,'t'))% check if excitation current has been defined
    flag = false;
    if ~isfield(D,'fcf')||~isfield(D,'Icf')||~isfield(D,'phf')
        error([ 'Please, define arrays of time waveform Ic and t,',...
                ' or frequencies, phases and current harmonics Icf, fcf and phf.']);
    end
else
    flag = true;
end

if flag  % current waveform in time has been provided
    % ensure column vectors
    D.t = D.t(:);
    D.Ic= D.Ic(:);
    
    % check vector sizes
    D.ni = numel(D.Ic);     % # of points in one cycle of current waveform
    if numel(D.t)~=D.ni
        error('Ensure that the number of elements in Ic is equal to those in t.');
    end
    
    % frequencies to evaluate
    nfDefined = false;
    if ~isfield(D,'nf')     % if # of freq to evaluate was not defined
        D.nf = 10;          % set default # of freq. to evaluate (D.nf)
    else                    % if D.nf was defined, but equals to 0
        if D.nf==0
            D.nf = 10;      % set default # of freq. to evaluate (D.nf)
        elseif D.nf>D.ni    % # of freq should not be larger than sample points
            error('Provide more sample points in Ic than freq. to evaluate.');
        else
            nfDefined = true;
        end
    end
    
    % determine harmonic spectrum -----------------------------------------
        
    % determine sampling period
    dt = (D.t(end)-D.t(1))/D.ni;    % sampling period [s]
    fs  = 1/dt;                     % sampling frequency [Hz]
    
    % compute FFT
    Y = fft(D.Ic,D.ni);             % compute FFT
    fsp = fs*(0:floor(D.ni/2))'/D.ni;% frequency spectrum [Hz]
    P2 = abs(Y/D.ni);               % calculate double sided spectrum (amplitude)
    ph2 = angle(Y);                 % phase information of the FFT
    P1 = P2(1:floor(D.ni/2)+1);     % calculate single sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);    % correct single sided spectrum
    ph1 = ph2(1:floor(D.ni/2)+1);   % single sided phase information of the FFT
    
    % get the largest harmonic currents to evaluate the model
    [Icsorted,idxs] = sort(P1,'descend');
    Icf = Icsorted(1:round(D.nf));
    % set each harmonic current as a fraction of the largest harmonic
    Ic_ratio = Icf/Icf(1);          % ratios of the largest current harmonic
    if ~nfDefined
        % evaluate only harmonics larger than 1% of the largest current harmonic
        if  length(Ic_ratio(Ic_ratio>0.01))<D.nf
            D.nf = length(Ic_ratio(Ic_ratio>0.01));
        end
    end
    D.Icf = Icsorted(1:round(D.nf));    % largest harmonic currents [A]
    D.fcf = fsp(idxs(1:round(D.nf)));   % freq. of largest harmonic currents [Hz]
    D.phf = ph1(idxs(1:round(D.nf)));   % phases of largest harmonic currents [rad]
else    % frequencies, phases and harmonic currents have been provided
    D.nf = length(D.Icf);% # of harmonic frequencies to evaluate
    % ensure column vectors
    D.Icf = D.Icf(:);   % current amplitude of each harmonic 
    D.fcf = D.fcf(:);   % frequency of each harmonic
    D.phf = D.phf(:);   % phase of each harmonic
    % create a virtual time vector, when this is not provided
    if ~isfield(D,'t')
        fnon0 = D.fcf(D.fcf~=0);
        fmin_non0 = min(fnon0);
        if isempty(fmin_non0)
            D.t   = transpose(linspace(0,1));
        else
            D.t   = transpose(linspace(0,1/fmin_non0));
        end
    end
end

%%
% compute harmonic current RMS value and radian frequencies
D.wn = 2*pi*D.fcf;                  % harmonic radian frequencies [rad/s]
D.irmsn = D.Icf;                    % get rms value of each harmonic current [A]
D.irmsn(D.fcf~=0) = D.irmsn(D.fcf~=0)/sqrt(2);

% compute time domain fitted current
mem = memory;
if numel(D.t)*D.nf>(mem.MaxPossibleArrayBytes/8)/2
    % compute fdfit by accumulation to save memory
    ifit = D.Icf(1)*cos(D.wn(1)*D.t+D.phf(1));
    for k = 2:nf
        ifit = ifit+D.Icf(k)*cos(D.wnf(k)*D.t+D.phf(k));
    end
else
    ifit = sum(D.Icf.*cos(D.wn*transpose(D.t)+D.phf),1); % fitted waveform
end

% if D.Ic was not provided, assign ifit to D.Ic
if ~isfield(D,'Ic'), D.Ic = ifit; end

% determine peak current
D.Ipk = max(D.Ic);                      % peak current [A]

%%
% assign output
varargout = {D};

%%
% plotting if figure number is parsed
if nargin>1
    % plot current excitation spectrum
    set(figure(fn),'color','white'); fn=fn+1;
    stem(D.fcf,D.Icf,'-');
    grid on; grid minor;
    xlabel('$f$ [Hz]','interpreter','latex');
    ylabel('$i_c(f)$ [A]','interpreter','latex');
    title('Current excitation frequency spectrum $i_c(f)$',...
        'interpreter','latex');
    
    % plot current excitation
    set(figure(fn),'color','white'); fn=fn+1;
    if flag
        plot(D.t,D.Ic,'-',D.t,ifit,'--');
        legend({'Input','Fit'},'interpreter','latex');
    else
        plot(D.t,ifit,'-');
    end
    grid on; grid minor;
    xlabel('$t$ [s]','interpreter','latex');
    ylabel('$i_c(t)$ [A]','interpreter','latex');
    title('Current excitation waveform $i_c(t)$ in one cycle',...
        'interpreter','latex');
    
    % assign output
    varargout = {D,fn};
end
end