clear all, clc, close all

%Initiation
%
%Constants of levitation medium
c0water=1510;         %Speed of sound in water     [m/s]1510
rho0water=997;        %Density of water            [kg/m^3]997

%Transducer array
f=1000000;            %Transducer frequency         [Hz]
a=0.0065;             %Piston radius                [m]
Radius=0.027;         %Radius of the cap and alignment height [m] 25mm for 43 deg
NR(1:4)=[10,0,0,0];    %Number of transducers i'th ring 
Rot(1:4)=[0,1,0,0];   %Ring i shited (1=yes; 0=no)
Alpha=75;             %Angle distribution w.r.t cap radius [deg]  (Default = 43-60(8)-75(10))
Trap=2;               %Choose "1" for vortex trap or "2" for twin trap

M=1;                  %Periodicity of the vortex                  (Default = 1)
% deg1=rad2deg(0);               %Vortex trap phases
% deg2=rad2deg(pi/5);
% deg3=rad2deg(2*pi/5);
% deg4=rad2deg(3*pi/5);
% deg5=rad2deg(4*pi/5);
% deg6=rad2deg(5*pi/5);
% deg7=rad2deg(6*pi/5);
% deg8=rad2deg(7*pi/5);
% deg9=rad2deg(8*pi/5);
% deg10=rad2deg(9*pi/5);


deg1=rad2deg(pi);               %Twin trap phases
deg2=rad2deg(pi);
deg3=rad2deg(pi);
deg4=rad2deg(pi);
deg5=rad2deg(pi);
deg6=rad2deg(0);
deg7=rad2deg(0);
deg8=rad2deg(0);
deg9=rad2deg(0);
deg10=rad2deg(0);

%Pressure constant
Measured_pressure=65000;        %[Pa]        1MHz measurement Imasonic
Vpp=25;                         %[Vpp]
Measured_distance=0.028;         %[m]
Vpp_actual=20;

%Define grid
gridsize=800;        %
Xlength=0.01;        %Half of the width of the plane [m]0.01
Zstart=0.015;        %Begin of imaging on z-axis     [m]0.015
dx=2*Xlength/gridsize;

%Ultraino spherical cap (x,y,z) location and phases of the Nt transducers
%
Anglespace=(Alpha*2*pi)/360;
ringangle(1:4)=[1:4]*Anglespace; %Angle of i'th ring w.r.t. array center [rad]
for l=1:4
Zi(l)=Radius-Radius*cos(ringangle(l));            %z-height of the i'th transducer ring   [m]         
end   
for l=1:4
AR(l)=Radius*sin(ringangle(l));                  %xy-plane radius of the i'th ring       [m]
end
[xarray,yarray,zarray,phi,Nt]=transducerloc(NR,AR,Zi,Rot,M,Trap,deg1,deg2,deg3,deg4,deg5,deg6,deg7,deg8,deg9,deg10);

%Acoustic wave complex pressure field calculation
%
lambda_water=c0water/f;                            %Wavelength in water            [m]
Theoretical_Focus=((2*a)^2)/(4*lambda_water)
k_water=(2*pi)/lambda_water;                       %Wave number for 40kHz in water [1/m]
Pressure_constant=Measured_pressure*Measured_distance/Vpp;
vn=(2*Pressure_constant*Vpp)/(i*k_water*rho0water*c0water*(a^2));      %Piston velocity [m/s]
Xloc=linspace(-Xlength,Xlength,gridsize);          %Define grid for x-axis
Yloc=linspace(-Xlength,Xlength,gridsize);          %Define grid for y-axis 
Zloc=linspace(Zstart,Zstart+2*Xlength,gridsize);   %Define grid for z-axis
Yaxis=0; Xaxis=0;
Pcomp_field_xz=zeros(gridsize,gridsize);           %xz-plane grid for complex pressure
Pcomp_field_xy=zeros(gridsize,gridsize);           %xy-plane grid for complex pressure
Pcomp_field_yz=zeros(gridsize,gridsize);           %yz-plane grid for complex pressure

%Calculation
for n=2:Nt                                         %Start with 2nd transducer, because origin no amplitude in vortex trap
for z=1:length(Zloc)
for x=1:length(Xloc)
   [d,theta]=distance_and_angle(Xloc(x),Yaxis,Zloc(z),xarray(n),yarray(n),zarray(n),Radius);
   Df_loc=directfunc(k_water,a,theta);
   Pcomp(x,z)=complex_pressure(k_water,rho0water,c0water,a,vn,Df_loc,d,phi(n));
   Pcomp_field_xz(x,z)=Pcomp_field_xz(x,z)+Pcomp(x,z);
end
end
end

%Extracting absolute pressure
Pcomp_field_xz=Pcomp_field_xz';
for x=1:gridsize
    for z=1:gridsize
        Ptot_xz(x,z)=abs(Pcomp_field_xz(x,z));
    end
end
[Ma,I]=max((Ptot_xz'));
[Pmax,Imax]=max(Ma)
Zplane=Zstart+(Imax*2*Xlength/gridsize);

%yz plane pressure
for n=2:Nt                                         %Start with 2nd transducer, because origin no amplitude in vortex trap
for z=1:length(Zloc)
for y=1:length(Yloc)
   %[d,theta]=distance_and_angle_yz(Xloc(x),Yloc(y),Zplane,xarray(n),yarray(n),zarray(n),Radius);
   [d,theta]=distance_and_angle_yz(Xaxis,Yloc(y),Zloc(z),xarray(n),yarray(n),zarray(n),Radius);
   Df_loc=directfunc(k_water,a,theta);
   Pcomp(y,z)=complex_pressure(k_water,rho0water,c0water,a,vn,Df_loc,d,phi(n));
   Pcomp_field_yz(y,z)=Pcomp_field_yz(y,z)+Pcomp(y,z);
end
end
end

%Extracting absolute pressure
Pcomp_field_yz=Pcomp_field_yz';
for y=1:gridsize
    for z=1:gridsize
        Ptot_yz(y,z)=abs(Pcomp_field_yz(y,z));
    end
end
[Ma2,I2]=max((Ptot_yz'));
[Pmax2,Imax2]=max(Ma2)
%Zplane=Zstart+(Imax*2*Xlength/gridsize);

%xy-plane pressure
for n=2:Nt                                         %Start with 2nd transducer, because origin no amplitude in vortex trap
for y=1:length(Yloc)
for x=1:length(Xloc)
   [d,theta]=distance_and_angle_xy(Xloc(x),Yloc(y),Zplane,xarray(n),yarray(n),zarray(n),Radius);
   Df_loc=directfunc(k_water,a,theta);
   Pcomp_xy(x,y)=complex_pressure(k_water,rho0water,c0water,a,vn,Df_loc,d,phi(n));
   Pcomp_field_xy(x,y)=Pcomp_field_xy(x,y)+Pcomp_xy(x,y);
end
end
end
Pcomp_field_xy=Pcomp_field_xy';
%Extracting absolute pressure
for x=1:gridsize
    for y=1:gridsize
        Ptot_xy(x,y)=abs(Pcomp_field_xy(x,y));
    end
end
[Ma3,I3]=max((Ptot_xy'));
[Pmax3,Imax3]=max(Ma3)
%%
%Visualising absolute pressure
figure(2),clf(2),
subplot(2,2,1), hold on
imagesc(1000*Xloc,1000*Zloc,0.001*Ptot_xz)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
colormap(hot)
colorbar
axis equal
%axis([-5,5,20,30])
axis([-10,10,15,35])
%title('xz-plane [kPa]')
%xlabel('X [mm]')
ylabel('Z [mm]')
hold off
subplot(2,2,3), hold on
imagesc(1000*Xloc,1000*Yloc,0.001*Ptot_xy)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
colormap(hot)
colorbar
axis equal
axis([-10,10,-10,10])
%axis([-5,5,-5,5])
%title('xy-plane [kPa]')
xlabel('X [mm]')
ylabel('Y [mm]')
hold off
%%
%Acoustic Radiation Force calculation
%
c0Glass=2350;         %Speed of sound in PS      [m/s]900 EPS 4540 Glass 2350 PS
rho0PS=1050;       %Density of PS             [kg/m^3]29 EPS 2500 Glass 1050 PS
dp=0.5;               %Particle diameter            [mm]

partdiam=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
for l=1:10
    for n=1:gridsize
Fg(l,n)=9.81*rho0PS*(4/3)*pi*(((partdiam(l))/2000)^3);
    end
end
%Potential wave energy field calculation
%Constants of particle

rp=dp/2000;                                        %Radius of the particle         [m] 
Uxz=potential_field(gridsize,Xlength,f,rp,rho0water,c0water,rho0PS,c0Glass,Ptot_xz,Pcomp_field_xz);
Uxy=potential_field(gridsize,Xlength,f,rp,rho0water,c0water,rho0PS,c0Glass,Ptot_xy,Pcomp_field_xy);
Uyz=potential_field(gridsize,Xlength,f,rp,rho0water,c0water,rho0PS,c0Glass,Ptot_yz,Pcomp_field_yz);
figure(5),clf(5), hold on
%subplot(1,2,1), hold on
imagesc(1000*Xloc,1000*Zloc,Uxz)
colormap(jet)
colorbar
axis equal
axis(1000*[-Xlength,Xlength,Zstart,Zstart+2*Xlength])
title('U in xz-plane [Nm]')
xlabel('x [mm]')
ylabel('z [mm]')
hold off
%subplot(1,2,2), hold on
%imagesc(1000*Xloc,1000*Yloc,Uxy)
%colormap(jet)
%colorbar
%axis equal
%axis(1000*[-Xlength,Xlength,-Xlength,Xlength])
%title('Acoustic potential wave field, U [Nm]')
%xlabel('x [mm]')
%ylabel('y [mm]')
%hold off

%Force vector field plot
%Phase
for x=1:gridsize
    for z=1:gridsize
        Pphase(x,z)=angle(Pcomp_field_xz(x,z));
    end
end
figure(4),clf(4), hold on
imagesc(1000*Xloc,1000*Zloc,Pphase)
colormap(hsv)
colorbar
axis equal
title('Acoustic pressure field phase [rad]')
xlabel('X [mm]')
ylabel('Z [mm]')
hold off

Zoomfactor=20;

Xvec=linspace(-1000*Xlength,1000*Xlength,(gridsize/Zoomfactor));
Yvec=linspace(-1000*Xlength,1000*Xlength,(gridsize/Zoomfactor));
Zvec=linspace(1000*Zstart,1000*Zstart+2*1000*Xlength,(gridsize/Zoomfactor));
% Fxz=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
% Fyz=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
% Fzz=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
% Fxy=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
% Fyy=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
% for x=3:((gridsize/Zoomfactor)-2)
%     for z=3:((gridsize/Zoomfactor)-2)
%         Fxz(z,x)=-(Uxz(Zoomfactor*z,Zoomfactor*x+Zoomfactor)-Uxz(Zoomfactor*z,Zoomfactor*x-Zoomfactor))/(2*dx*Zoomfactor);
%         Fzz(z,x)=-(Uxz(Zoomfactor*z+Zoomfactor,Zoomfactor*x)-Uxz(Zoomfactor*z-Zoomfactor,Zoomfactor*x))/(2*dx*Zoomfactor);
%     end
% end
% for y=3:((gridsize/Zoomfactor)-2)
%     for z=3:((gridsize/Zoomfactor)-2)
%         Fyz(y,z)=-(Uyz(Zoomfactor*z,Zoomfactor*y+Zoomfactor)-Uyz(Zoomfactor*z,Zoomfactor*y-Zoomfactor))/(2*dx*Zoomfactor);
%         Fzz(z,y)=-(Uyz(Zoomfactor*z+Zoomfactor,Zoomfactor*y)-Uxz(Zoomfactor*z-Zoomfactor,Zoomfactor*y))/(2*dx*Zoomfactor);
%     end
% end
% for x=3:((gridsize/Zoomfactor)-2)
%     for y=3:((gridsize/Zoomfactor)-2)
%         Fxy(y,x)=-(Uxy(Zoomfactor*y,Zoomfactor*x+Zoomfactor)-Uxy(Zoomfactor*y,Zoomfactor*x-Zoomfactor))/(2*dx*Zoomfactor);
%         Fyy(y,x)=-(Uxy(Zoomfactor*y+Zoomfactor,Zoomfactor*x)-Uxy(Zoomfactor*y-Zoomfactor,Zoomfactor*x))/(2*dx*Zoomfactor);
%     end
% end
%figure(6), clf(6)
%subplot(1,2,1), hold on
%quiver(Xvec,Zvec,Fxz,Fzz,'k');
%axis equal
%title('Acoustic radiation force vector field F_j [N]')
%xlabel('x [mm]')
%ylabel('z [mm]')
%hold off
%subplot(1,2,2), hold on
%quiver(Xvec,Zvec,Fxy,Fyy,'k');
%axis equal
%title('Acoustic radiation force vector field F_j [N]')
%xlabel('x [mm]')
%ylabel('y [mm]')
%hold off


%Trapping force plots
Fxtrap=0*[1:gridsize];
Fytrap=0*[1:gridsize];
Fztrap=0*[1:gridsize];
Xtrap=linspace(-1000*Xlength,1000*Xlength,(gridsize));
Ytrap=linspace(-1000*Xlength,1000*Xlength,(gridsize));
Ztrap=linspace(0,2*1000*Xlength,(gridsize));

for x=5:(gridsize-5)
        Fxtrap(x)=1000000*-(Uxz(Imax,x+1)-Uxz(Imax,x-1))/(2*dx);
        %Fxtrap(x)=1000000*-(Uxz((gridsize/2),x+1)-Uxz((gridsize/2),x-1))/(2*dx);
end
for y=5:(gridsize-5)
        %Fytrap(y)=-1000000*(Uyz(z+1,(gridsize/2))-Uyz(z-1,(gridsize/2)))/(2*dx);
        Fytrap(y)=1000000*(Uyz(Imax2,y+1)-Uyz(Imax2,y-1))/(2*dx);
end
for z=5:(gridsize-5)
        Fztrap(z)=-1000000*(Uxz(z+1,(0*100+gridsize/2))-Uxz(z-1,(0*100+gridsize/2)))/(2*dx);
end

Fgrav=Fg(1,:);

figure(2)
%figure(7), clf(7), hold on
subplot(2,2,4), hold on
plot(Xtrap,Fxtrap,'-r','LineWidth',2)
hold on
plot([-0.337,-0.337],[-10,10],'--k','linewidth',2) %68.98 [-0.413,-0.413]
hold on
plot([0.337,0.337],[-10,10],'--k','linewidth',2) %68.98
hold on
% plot(Xtrap,Fgrav,'-k')
% hold on
% plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
xlabel('X [mm]')
axis([-4,4,-1.5,1.5])
grid on
ylabel('F_{rad,lat} [\muN]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
hold off
%Fgrav(1)

%figure(8), clf(8), hold on
subplot(2,2,2), hold on
plot(2*Ztrap,Fztrap,'-r','LineWidth',2)             % Multiplying with 2 for accuracy
hold on
% plot([-0.337,-0.337],[-10,10],'--k','linewidth',2) %68.98 [-0.413,-0.413]
% hold on
% plot([0.337,0.337],[-10,10],'--k','linewidth',2) %68.98
% hold on
% plot(Xtrap,Fgrav,'-k')
% hold on
% plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
xlabel('Z [mm]')
%axis([-4, 4, -0.02, 0.02])%[-4,4,-5,5])
grid on
ylabel('F_{rad,axial} [\muN]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
hold off

%% Plotting all components of forces together
figure(3)
%figure(7), clf(7), hold on
subplot(2,1,2), hold on
plot(Xtrap,Fxtrap,'<r','LineWidth',1)
hold on
plot([-0.337,-0.337],[-10,10],'--k','linewidth',2) %68.98 [-0.413,-0.413]
hold on
plot([0.337,0.337],[-10,10],'--k','linewidth',2) %68.98
hold on
% plot(Xtrap,Fgrav,'-k')
% hold on
% plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
xlabel('X [mm]')
axis([-4,4,-1.5,1.5])
grid on
ylabel('F_{rad,lat} [\muN]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 16)
hold off
%Fgrav(1)

%figure(8), clf(8), hold on
% subplot(3,1,2), hold on
% plot(Ytrap,Fytrap,'-r','LineWidth',2)
% hold on
% plot([-0.337,-0.337],[-10,10],'--k','linewidth',2) %68.98 [-0.413,-0.413]
% hold on
% plot([0.337,0.337],[-10,10],'--k','linewidth',2) %68.98
% hold on
% plot(Xtrap,Fgrav,'-k')
% hold on
% plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
% xlabel('x [mm]')
% %axis([-4, 4, -0.02, 0.02])%[-4,4,-5,5])
% grid on
% ylabel('F_{rad,y} [\muN]')
% set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
% hold off

%figure(8), clf(8), hold on
subplot(2,1,1), hold on
plot(Zstart*ones(length(Ztrap),1)'+2*Ztrap,1000*Fztrap,'>r','LineWidth',1)  % multiplying a scalar 2 to Z axis check with RJ
hold on
% plot([-0.337,-0.337],[-10,10],'--k','linewidth',2) %68.98 [-0.413,-0.413]
% hold on
% plot([0.337,0.337],[-10,10],'--k','linewidth',2) %68.98
% hold on
% plot(Xtrap,Fgrav,'-k')
% hold on
% plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
xlabel('Z [mm]')
axis([0, 40, -5, 5])%[-4,4,-5,5])
grid on
ylabel('F_{rad,axial} [nN]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 16)
hold off


%figure(8), clf(8), hold on
%plot(Fztrap,Ztrap,'-k')
%plot(Fg(3,:),Ztrap,'-r')
   %plot(0*[1:length(Xtrap)],Ztrap,'-k')
%title('Axial trapping force at trapping center')
%xlabel('Trapping force F_z [N]')
%ylabel('z-axis [mm]')
%hold off
%Fmax=max(Fxtrap)
min_R=sin((pi-Anglespace)/2)*0.017/(sin(Anglespace));
min_AR=0.017/(2*sin(pi/NR(1)))
Actual_AR_1=AR(1)
%min_AR_16=sin((pi-(2*pi/16))/2)*0.025/(sin((2*pi/16)))
%min_AR_21=sin((pi-(2*pi/21))/2)*0.025/(sin((2*pi/21)))
%min_AR_28=sin((pi-(2*pi/28))/2)*0.025/(sin((2*pi/28)))
Fz_A = max(Fztrap)
Fx_A = max(Fxtrap)