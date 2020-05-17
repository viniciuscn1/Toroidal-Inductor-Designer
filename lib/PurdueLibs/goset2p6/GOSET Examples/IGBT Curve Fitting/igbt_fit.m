% this rountine returns the fitness of the IGBT conduction
% loss parameters based on the mean percentage error

function fitness = igbt(parameters,data,fignum)

% assign genes to parameters
a = parameters(1);
b = parameters(2);
c = parameters(3);

vpred = a*data.i+(b*data.i).^c;
error = abs(1-vpred./data.v);
fitness = 1.0/(1.0e-6+mean(error));

if nargin>2

   figure(fignum);
   Npoints=200;
   ip=linspace(0,max(data.i),Npoints);
   vp=a*ip+(b*ip).^c;
   plot(data.i,data.v,'bx',ip,vp,'r')
   title('Voltage Versus Currrent');
   xlabel('Current, A');
   ylabel('Voltage, V');
   legend({'Measured','Fit'}); 

   figure(fignum+1);
   Npoints=200;
   plot(data.i,data.v.*data.i,'bx',ip,vp.*ip,'r')
   title('Power Versus Currrent');
   xlabel('Current, A');
   ylabel('Power, V');
   legend({'Measured','Fit'}); 

end