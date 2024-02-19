function [ dNxi, dNeta ] = func_dNilin(xi,eta )


dNxi  = zeros(3,1);
dNeta = zeros(3,1);


dNxi(1) = -1;
dNxi(2) = 1;
dNxi(3) = 0;

dNeta(1) = -1;
dNeta(2) = 0;
dNeta(3) = 1;


end