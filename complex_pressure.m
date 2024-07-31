function [Pcomp]=complex_pressure(k,rho,c0,a,vn,Df_loc,d,phi)
Pcomp=(1/2)*i*k*rho*c0*(a^2)*vn*Df_loc*(1/d)*exp(i*(phi+k*d));
end
