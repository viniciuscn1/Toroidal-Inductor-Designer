function bodeplot(fignum,f,varargin)
% BODEPLOT plots the magnitude of a quantity in dB versus f on a semilog axis and
% the angle of a quantity in degrees versus f on a semilog axis.
%
% bodeplot(fignum,f,x1,s1,[x2,s2],[x3,s3],...)
%
% Inputs:
% fignum  = figure number
% f       = frequency, Hz
% x1      = data curve 1
% s1      = non-optional plot descriptor for data curve 1
% [x2,s2] = optional second curve and descriptor
% [x3,s3] = optional third curve and descriptor
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% check inputf
if nargin<4
   error('Error in magplot: insufficient number of arguments');
end
if (mod((nargin),2)~=0)
   error('Error in magplot: need odd number of arguments');
end

% compute number of curves
Nc=round((nargin-2)/2);

% set up figure
figure(fignum);
h1=subplot(2,1,1);

% draw first curve
x_dB=20.0*log10(abs(varargin{1}));
s=varargin{2};
semilogx(f,x_dB,s);
% draw remaining curves
if Nc>1
   hold on;
   for  i=2:Nc,
      x_dB=20.0*log10(abs(varargin{2*i-1}));
      s=varargin{2*i};
      semilogx(f,x_dB,s);
  end
  hold off;
end
grid on;
ylabel('magnitude, dB');


h2=subplot(2,1,2);
x_deg=angle(varargin{1})*180/pi;
s=varargin{2};
semilogx(f,x_deg,s);
grid on;
ylabel('phase, degrees');
xlabel('frequency, Hz');
% draw remaining curves
if Nc>1
   hold on;
   for  i=2:Nc,
      x_deg=angle(varargin{2*i-1})*180/pi;
      s=varargin{2*i};
      semilogx(f,x_deg,s);
   end
   hold off;
end

% label curves
subplot(h1);