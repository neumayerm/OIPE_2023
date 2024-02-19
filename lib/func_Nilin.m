function [ N ] = func_Nilin(xi,eta )


N = zeros(3,1);

 N(1) = 1 - xi - eta;
 N(2) = xi;
 N(3) = eta;

end

