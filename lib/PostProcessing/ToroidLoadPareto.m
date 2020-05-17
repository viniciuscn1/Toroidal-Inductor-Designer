function sol = ToroidLoadPareto(sol,filename,desc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load optimization results file for post-processing.
%                                                                         %
% Written by Vinicius Nascimento                                          %
% viniciuscn1@gmail.com
%                                                                         %
% Current revision date: 01-10-2020                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load file
data = load(filename);

% save design specs structure
sol(1).n            = sol(1).n+1;
sol(sol(1).n).D     = data.D;
sol(sol(1).n).bi    = data.bi;
sol(sol(1).n).desc  = desc;
sol(sol(1).n).GAP   = data.GAP;
sol(sol(1).n).GAS   = data.GAS;
sol(sol(1).n).fP    = data.fP;
sol(sol(1).n).GAP.rp_lvl = 1;

% extract information from parameters
sol(sol(1).n).NS = size(sol(sol(1).n).bi,2);
sol(sol(1).n).res = ToroidFit(sol(sol(1).n).bi(:,1),sol(sol(1).n).D,0);
for k = 1:sol(sol(1).n).NS
   f = ToroidFit(sol(sol(1).n).bi(:,k),sol(sol(1).n).D,0);
   sol(sol(1).n).res(k) = f;
   sol(sol(1).n).rc(k,:) = sol(sol(1).n).res(k).rperf;
   sol(sol(1).n).perf(k,:) = sol(sol(1).n).res(k).perf;
   sol(sol(1).n).lc(k)  = sol(sol(1).n).res(k).par.lc;
   sol(sol(1).n).ra(k) = sol(sol(1).n).res(k).par.ra;
   sol(sol(1).n).rci(k) = sol(sol(1).n).res(k).par.rci;
   sol(sol(1).n).rco(k) = sol(sol(1).n).res(k).par.rco;
   sol(sol(1).n).N(k) = sol(sol(1).n).res(k).W.N;
   sol(sol(1).n).awg(k) = str2double(sol(sol(1).n).res(k).W.desc);
   sol(sol(1).n).fmu(k,:) = transpose(sol(sol(1).n).res(k).par.fmu);
end

% sort design information by mass------------------------------------------
[sol(sol(1).n).mass,sol(sol(1).n).idx] = ...
    sort(transpose([sol(sol(1).n).res(:).ML]),'ascend');
sol(sol(1).n).loss = transpose([sol(sol(1).n).res(sol(sol(1).n).idx).PL]);
sol(sol(1).n).vol = transpose([sol(sol(1).n).res(sol(sol(1).n).idx).VL]);
sol(sol(1).n).lwr  = transpose([sol(sol(1).n).res(sol(sol(1).n).idx).lwr]);
sol(sol(1).n).aL   = transpose([sol(sol(1).n).res(sol(sol(1).n).idx).aL]);
sol(sol(1).n).Twmx = transpose([sol(sol(1).n).res(sol(sol(1).n).idx).Twmx]);
sol(sol(1).n).bi   = sol(sol(1).n).bi(:,sol(sol(1).n).idx);
sol(sol(1).n).id   = 1:length(sol(sol(1).n).idx);
sol(sol(1).n).rc   = sol(sol(1).n).rc(sol(sol(1).n).idx,:);
sol(sol(1).n).perf = sol(sol(1).n).perf(sol(sol(1).n).idx,:);

end