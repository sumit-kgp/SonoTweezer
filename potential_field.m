function [U]=potential_field(gridsize,Xlength,f,rp,rho0,c0,rho0part,c0part,Ptot,Pcomp_field)
U=zeros(gridsize,gridsize);             %Define xz-plane grid for potential energy
%constants
alpha=(1/3)*pi*(rp^3)*((1/(rho0*c0^2))-(1/(rho0part*c0part^2)));
beta=pi*(rp^3)*((-rho0+rho0part)/(((2*pi*f)^2)*rho0*(rho0+2*rho0part)));
dx=2*Xlength/gridsize;
%Calculation
for x=3:(gridsize-2)
    for z=3:(gridsize-2)
        dpdx=(1/(12*dx))*(-Pcomp_field(z,x+2)+8*Pcomp_field(z,x+1)-8*Pcomp_field(z,x-1)+Pcomp_field(z,x-2));
        dpdz=(1/(12*dx))*(-Pcomp_field(z+2,x)+8*Pcomp_field(z+1,x)-8*Pcomp_field(z-1,x)+Pcomp_field(z-2,x));
        %dpdx=(1/(2*dx))*(Pcomp_field(z,x+1)-Pcomp_field(z,x-1));
        %dpdz=(1/(2*dx))*(Pcomp_field(z+1,x)-Pcomp_field(z-1,x));
U(z,x)=alpha*(Ptot(z,x)^2)-beta*((abs(dpdz))^2+(abs(dpdx))^2);
    end
end
% alpha
% beta
% rt =alpha/beta
end
