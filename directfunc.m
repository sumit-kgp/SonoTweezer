function [Df] = directfunc(k,a,theta)
Df=(2*besselj(1,(k*a*sin(theta))))/(k*a*sin(theta));
end