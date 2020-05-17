%% ToroidPostProcess.m
% This script loads a result optimization file, and plot the Pareto Optimal
% it allows for the selection a specific design for visualization and
% analysis.
% 
% Written: Vinicius C. do Nascimento
%          viniciuscn1@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear workspace
close all;
clear;
clc;

% add libs
addpath(genpath([pwd,'\lib']));

% define local variables
fn = 1;                     % figure # control variable

% define local constants
um = 1e-6;                  % convert micrometer to meter [m]
mm = 1e-3;                  % convert milimeter to meter [m]

%%
% load results files
sol.n = 0;
sol = ToroidLoadPareto(sol,'Ip_15A_1000Hz_130oC_f3_results_Tamb40_nspc1','test');

% plot the GA results
distplot(100,sol(1).fP,1,sol(1).GAP);
reportplot(sol(1).GAP,sol(1).GAS,sol(1).fP);

%%
% look at a particular solution--------------------------------------------
mass_target = 0.2255;       % target mass in the Pareto optimal front
[~,sol(1).didx] = min(abs(sol(1).mass-mass_target));

if (sol(1).didx>0)
   if (sol(1).didx>length(sol(1).id))
       disp('No such design');
   else
       % show location of design on pareto optimal front
       set(figure(fn),'color','white'); fn=fn+1;
       if sol(1).D.nobj==3
           plot3(sol(1).mass,sol(1).loss,sol(1).vol,'x-',...
                sol(1).mass(sol(1).didx),sol(1).loss(sol(1).didx),...
                sol(1).vol(sol(1).didx),'ro');
            zlabel('Volume [$m^3$]','interpreter','latex');
       elseif sol(1).D.nobj==2
           plot(sol(1).mass,sol(1).loss,'x-',...
                sol(1).mass(sol(1).didx),sol(1).loss(sol(1).didx),'ro');
       end
       xlabel('Mass [kg]','interpreter','latex');
       ylabel('Loss [W]','interpreter','latex');
       grid on; grid minor;
       title('Pareto-Optimal Front','interpreter','latex');
       legend({'All Solutions',['Design ' num2str(sol(1).didx)]},...
           'interpreter','latex','location','best');
       % max temperature vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).Twmx-273.15,'*',...
           [sol(1).mass(1),sol(1).mass(end)],...
           ones(2,1)*sol(1).D.Twmxa-273.15,'--');
       title('Max. Winding Temperature vs Mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$T_{wmx}$, $^{\circ}C$','interpreter','latex');
       grid on; grid minor;
       ylim([0,sol(1).D.Twmxa-273.15]*1.2);
       % aspect ratio vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).aL,'-',...
           [sol(1).mass(1),sol(1).mass(end)],...
           ones(2,1)*sol(1).D.amxa,'--');
       title('Aspect Ratio vs Mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$a_{L}$','interpreter','latex');
       grid on; grid minor;
       % total wire length vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).lwr,'-',...
           [sol(1).mass(1),sol(1).mass(end)],...
           ones(2,1)*9.5,'--');
       title('Total wire length vs Mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$l_{wr}$','interpreter','latex');
       grid on; grid minor;
       % axial length vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).lc(sol(1).idx)/mm);
       title('Core Axial Length $l_c$ vs Mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$l_c$, $mm$','interpreter','latex');
       grid on; grid minor;
       % # of turns vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).N(sol(1).idx));
       title('Number of turns vs mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('N','interpreter','latex');
       grid on; grid minor;
       % AWG vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).awg(sol(1).idx));
       title('AWG vs mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('AWG','interpreter','latex');
       grid on; grid minor;
       % inner and outer radius vs mass
       set(figure(fn),'color','white'); fn=fn+1;
       plot(sol(1).mass,sol(1).ra(sol(1).idx)/mm,...
            sol(1).mass,sol(1).rci(sol(1).idx)/mm,...
            sol(1).mass,sol(1).rco(sol(1).idx)/mm);
       title('Core Inner and Outer Radii vs Mass','interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$r_c$, $mm$','interpreter','latex');
       legend({'$r_a$','$r_{ci}$','$r_{co}$'},'interpreter','latex','location','best');
       grid on; grid minor;
       % permeability profile for each design
       np = size(sol(1).perf,2);
       set(figure(fn),'color','white'); fn=fn+1;
       if length(sol(1).fmu(1,:))==1
           nmu = 1;
       else
           nmu = length(sol(1).fmu(1,:))-1; % number of permeability sections
       end
       x = zeros(sol(1).NS,np);
       y = zeros(sol(1).NS,np);
       z = zeros(sol(1).NS,np);
       for k = 1:sol(1).NS
           x(k,:) = sol(sol(1).n).mass(k)*ones(1,np);
           y(k,:) = sol(1).rc(k,:);
           z(k,:) = sol(1).perf(k,:)*sol.D.mp.mur;
       end              
       surf(x,y,z,'EdgeAlpha',0);
       shading interp;
       camlight right;
       lighting gouraud; 
       hcb = colorbar;
       title('Permeability function across radius for each design',...
           'interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$r$, $mm$','interpreter','latex');
       zlabel('$\mu_r$','interpreter','latex');
       grid on; grid minor;
       view(45,25);
       title(hcb,'$\mu_r$','interpreter','latex');
       
       % plot 2D version of permeability profiles
       set(figure(fn),'color','white'); fn=fn+1;
       plt = pcolor(x,y,z);
       set(plt,'FaceColor','interp','EdgeColor','interp');
       hold on;
       prci= plot(x(:,1),y(:,1),'k--','LineWidth',1);
       prco= plot(x(:,end),y(:,end),'k-','LineWidth',1);
       hold off;
       hcb = colorbar;
       title(hcb,'$\mu_r$','interpreter','latex');
       title('Permeability function across radius for each design',...
           'interpreter','latex');
       xlabel('Mass, kg','interpreter','latex');
       ylabel('$r$, $mm$','interpreter','latex');
       xlim([0.2,0.5]);
       ylim([0.01,0.035]);
       legend([prci,prco],{'$r_{ci}$','$r_{co}$'},...
           'interpreter','latex','location','northwest');
       grid on; grid minor;
       
       % give a detailed report on the design
       ToroidEval(sol(1).bi(:,sol(1).didx),sol(1).D,fn);
       
    end
end
