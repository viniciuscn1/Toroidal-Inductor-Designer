function c=gte(x,xmn)
% GTE        greater-than-or-equal to function
%
% c=gte(x,xmn)
%
% Inputs:
% x        = a quantitity 
% xmn      = minimum allowed value
%
% Outputs:
% c        = constraint variable - 1 if x>=xmn, 0<c<1 if x<xmn

   if (x>=xmn)
      c=1;
   else
      c=1/(1+xmn-x);
   end

end