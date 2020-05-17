function c=lte(x,xmx)
% LTE      less-than-or-equal to function
%
% c=lte(x,xmx)
%
% Inputs:
% x        = a quantitity 
% xmx      = maximum allowed value
%
% Outputs:
% c        = constraint variable - 1 if x<=xmx, 0<c<1 if x>xmx

   if (x<=xmx)
      c=1;
   else
      c=1/(1+x-xmx);
   end
   
end