% Power diode V-I curve fitting

% Voltage-current measuerment data for power diode in a Fuji 6MBI30L-060 
load igbt_data.mat

% plot the raw data
figure(1);
plot(igbt_data.i,igbt_data.v,'x');
xlabel('Current, A');
ylabel('Voltage, V');
title('Static IGBT Characteristics at 25 C');

% Initialize the parameters
GAP = gapdefault(1,1,500,100);
GAP.mg_nreg = 4;  
GAP.mg_tmig = 20;  
GAP.mg_pmig = 0.05;

% Define gene parameters 
%          min  max  chrm chrm     par  
%          val  val  type  id       #   description
GAP.gd = [ 1e-8 1e+0   3  1; ...  % 1  a
           1e-6 1e+3   3  1; ...  % 2  b
           1e-3 1e+0   3  1];     % 3  c

% Execute GOSET
[P,GAS,bestx] = gaoptimize(@igbt_fit,GAP,igbt_data);

% Plot the measured V-I data and the estimated V-I curve
igbt_fit(bestx,igbt_data,4);