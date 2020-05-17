%% ToroidComparePareto.m
% This script loads result optimization files, and plot Pareto Optimal
% curves together for comparison purposes.
% 
% Written: Vinicius C. do Nascimento
%          viniciuscn1@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% clear workspace
close all;
clear;
clc;

% load Pareto fronts
sol.n = 0;
sol = ToroidLoadPareto(sol,'Ip_15A_1000Hz_130oC_f1_results_Tamb40_nspc1','Constant$ \mu$');
sol = ToroidLoadPareto(sol,'Ip_15A_1000Hz_130oC_f2_results_Tamb40_nspc1','Linear $\mu$');
sol = ToroidLoadPareto(sol,'Ip_15A_1000Hz_130oC_f3_results_Tamb40_nspc1','Quadratic $\mu$');
sol = ToroidLoadPareto(sol,'Ip_15A_1000Hz_130oC_f4_results_Tamb40_nspc1','Cubic $\mu$');

%%
% show pareto optimal front
symbols = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
if length(sol)>strlength(symbols)
    n = length(sol)/strlength(symbols);
    for k = 1:ceil(log2(n))
        symbols = strcat(symbols,symbols);
    end
end
fn = 2;
set(figure(fn),'color','white'); fn=fn+1;
nsol = length(sol);
cap = cell(nsol,1);
for k = 1:nsol
    plot(sol(k).mass,[sol(k).res(sol(k).idx).PL],symbols(k));
    cap{k} = sol(k).desc;
    hold on;
end
hold off;
xlabel('Mass, kg','interpreter','latex');
ylabel('Loss, W','interpreter','latex');
title('Pareto-Optimal Front','interpreter','latex');
legend(cap,'interpreter','latex','location','best');
grid on; grid minor;

%%
% compare temperature
set(figure(fn),'color','white'); fn=fn+1;
nsol = length(sol);
cap = cell(nsol,1);
for k = 1:nsol
    plot(sol(k).mass,sol(k).Twmx-273.15,symbols(k));
    cap{k} = sol(k).desc;
    hold on;
end
hold off;
title('Max. Winding Temperature vs Mass','interpreter','latex');
xlabel('Mass, kg','interpreter','latex');
ylabel('$T_{wmx}$, $^{\circ}C$','interpreter','latex');
legend(cap,'interpreter','latex','location','best');
grid on; grid minor;

%%
% Pareto+Temperature
symbols = ['o','s','p','d','h','+','*','^','v','>','<','.','x'];
colors  = [1 0 0;
           0 0 0;
           0.6350 0.0780 0.1840;
           0 0 1;
           0.8500 0.3250 0.0980;
           0.4940 0.1840 0.5560;
           0.6350 0.0780 0.1840];
style   = {'-','-.','--',':'};
set(figure(fn),'color','white'); fn=fn+1;
nsol = length(sol);
cap = cell(nsol,1);
s = zeros(nsol,1);
p = zeros(nsol,1);
for k = 1:nsol
    s(k) = scatter(sol(k).mass,[sol(k).res(sol(k).idx).PL],...
        30,sol(k).Twmx-273.15,symbols(k),'filled');
    hold on;
    cap{k} = sol(k).desc;
    p(k) = plot(sol(k).mass,[sol(k).res(sol(k).idx).PL],...
        'LineStyle',style{k},'color',colors(k,:),'LineWidth',2);
end
xlabel('Mass, kg','interpreter','latex');
ylabel('Loss, W','interpreter','latex');
title('Pareto-Optimal Front','interpreter','latex');
grid on; grid minor;
hcb = colorbar;
title(hcb,'Peak Temperature [$^{\circ}$C]','interpreter','latex');
ylim([12,24]);
xlim([0.1,0.5]);
caxis([70,130]);
% selected design
pidx = 3;
sidx = 15;
ps = plot(sol(pidx).mass(sidx),sol(pidx).loss(sidx),'o','LineWidth',2);
txtmn = ['Design ',num2str(sidx),' \rightarrow '];
text(sol(pidx).mass(sidx),sol(pidx).loss(sidx),txtmn,'HorizontalAlignment','right');
legend(p,cap,'interpreter','latex','location','northwest');