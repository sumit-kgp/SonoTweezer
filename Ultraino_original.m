clear all, clc, close all

%Initiation
%
%Constants of levitation medium
c0air=343;            %Speed of sound in air        [m/s]
rho0air=1.18;         %Density of air               [kg/m^3]

%Transducer array
f=40000;            %Transducer frequency         [Hz]
off = 5;
a=0.0045;             %Piston radius                [m]
Radius=0.061;         %Radius of the cap and alignment height [m] 
NR(1:4)=[8,12,24,0];    %Number of transducers i'th ring 
Rot(1:4)=[4,6,12,0]+[off,off,off,off];   %Ring i shited (1=yes; 0=no) [1,0,1,0]; [6,8,14,0]+off-3 - all horizontal
Alpha=15;             %Angle distribution w.r.t cap radius [deg] 
Trap=1;               %Choose "1" for vortex trap or "2" for twin trap (SonoTweezer) or "3" for twin trap otherwise

M=1;                  %Periodicity of the vortex                  
deg1=180;               %Twin trap phases
deg2=180;
deg3=180;
deg4=180;
deg5=0;
deg6=0;
deg7=0;
deg8=0;

%Define grid
gridsize=800;        %
Xlength=0.04;        %Half of the width of the plane [m]
Zstart=0.012;        %Begin of imaging on z-axis     [m] z =0.022
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
[xarray,yarray,zarray,phi,Nt]=transducerloc(NR,AR,Zi,Rot,M,Trap,deg1,deg2,deg3,deg4,deg5,deg6,deg7,deg8);

% checking transducer distribution
figure(8)
plot3(xarray,yarray,zarray,'-o')
%%
%Acoustic wave complex pressure field calculation
%
lambda_air=c0air/f;                            %Wavelength            [m]
Theoretical_Focus=((2*a)^2)/(4*lambda_air)
k_air=(2*pi)/lambda_air;                       %Wave number  [1/m]
Pressure_constant=0.17;
Vpp=17;
vn=(2*Pressure_constant*Vpp)/(i*k_air*rho0air*c0air*(a^2));      %Piston velocity [m/s]
Xloc=linspace(-Xlength,Xlength,gridsize);          %Define grid for x-axis
Yloc=linspace(-Xlength,Xlength,gridsize);          %Define grid for x-axis 
Zloc=linspace(Zstart,Zstart+2*Xlength,gridsize);   %Define grid for z-axis
Yaxis=0;
Pcomp_field_xz=zeros(gridsize,gridsize);           %xz-plane grid for complex pressure
Pcomp_field_xy=zeros(gridsize,gridsize);           %xy-plane grid for complex pressure


%Calculation
for n=2:Nt                                         %Start with 2nd transducer, because origin no amplitude in vortex trap
for z=1:length(Zloc)
for x=1:length(Xloc)
   [d,theta]=distance_and_angle(Xloc(x),Yaxis,Zloc(z),xarray(n),yarray(n),zarray(n),Radius);
   Df_loc=directfunc(k_air,a,theta);
   Pcomp(x,z)=complex_pressure(k_air,rho0air,c0air,a,vn,Df_loc,d,phi(n));
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
%xy-plane pressure
for n=2:Nt                                         %Start with 2nd transducer, because origin no amplitude in vortex trap
for y=1:length(Yloc)
for x=1:length(Xloc)
   [d,theta]=distance_and_angle_xy(Xloc(x),Yloc(y),Zplane,xarray(n),yarray(n),zarray(n),Radius);
   Df_loc=directfunc(k_air,a,theta);
   Pcomp_xy(x,y)=complex_pressure(k_air,rho0air,c0air,a,vn,Df_loc,d,phi(n));
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
%%
%Visualising absolute pressure
figure(2),clf(2),
subplot(2,1,1), hold on
imagesc(1000*Xloc,1000*Zloc,Ptot_xz)
colormap(hot)
colorbar
axis equal
%title('xz-plane [Pa]')
axis([-10,10,15,35])
xlabel('X [mm]')
ylabel('Z [mm]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
hold off
subplot(2,1,2), hold on
imagesc(1000*Xloc,1000*Yloc,Ptot_xy)
colormap(hot)
colorbar
axis equal
%title('xy-plane [Pa]')
axis([-10,10,-10,10])
xlabel('X [mm]')
ylabel('Y [mm]')
set(gca, 'FontName', 'CMU Serif', 'FontSize', 24)
hold off
%%
%Acoustic Radiation Force calculation
%
c0EPS=900;         %Speed of sound in PS      [m/s]900 EPS 4540 Glass 2350 PS
rho0EPS=29;       %Density of PS             [kg/m^3]29 EPS 2500 Glass 1050 PS
dp=1;               %Particle diameter            [mm]

partdiam=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
for l=1:10
    for n=1:gridsize
Fg(l,n)=9.81*rho0EPS*(4/3)*pi*(((partdiam(l))/2000)^3);
    end
end
%Potential wave energy field calculation
%Constants of particle

rp=dp/2000;                                        %Radius of the particle         [m] 
Uxz=potential_field(gridsize,Xlength,f,rp,rho0air,c0air,rho0EPS,c0EPS,Ptot_xz,Pcomp_field_xz);
Uxy=potential_field(gridsize,Xlength,f,rp,rho0air,c0air,rho0EPS,c0EPS,Ptot_xy,Pcomp_field_xy);
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
xlabel('x [mm]')
ylabel('z [mm]')
hold off

Zoomfactor=20;

Xvec=linspace(-1000*Xlength,1000*Xlength,(gridsize/Zoomfactor));
Zvec=linspace(1000*Zstart,1000*Zstart+2*1000*Xlength,(gridsize/Zoomfactor));
Fxz=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
Fzz=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
Fxy=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
Fyy=zeros((gridsize/Zoomfactor),(gridsize/Zoomfactor));
for x=3:((gridsize/Zoomfactor)-2)
    for z=3:((gridsize/Zoomfactor)-2)
        Fxz(z,x)=-(Uxz(Zoomfactor*z,Zoomfactor*x+Zoomfactor)-Uxz(Zoomfactor*z,Zoomfactor*x-Zoomfactor))/(2*dx*Zoomfactor);
        Fzz(z,x)=-(Uxz(Zoomfactor*z+Zoomfactor,Zoomfactor*x)-Uxz(Zoomfactor*z-Zoomfactor,Zoomfactor*x))/(2*dx*Zoomfactor);
    end
end
for x=3:((gridsize/Zoomfactor)-2)
    for y=3:((gridsize/Zoomfactor)-2)
        Fxy(y,x)=-(Uxy(Zoomfactor*y,Zoomfactor*x+Zoomfactor)-Uxy(Zoomfactor*y,Zoomfactor*x-Zoomfactor))/(2*dx*Zoomfactor);
        Fyy(y,x)=-(Uxy(Zoomfactor*y+Zoomfactor,Zoomfactor*x)-Uxy(Zoomfactor*y-Zoomfactor,Zoomfactor*x))/(2*dx*Zoomfactor);
    end
end
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
Fztrap=0*[1:gridsize];
Xtrap=linspace(-1000*Xlength,1000*Xlength,(gridsize));
Ztrap=linspace(0,2*1000*Xlength,(gridsize));

for x=5:(gridsize-5)
        Fxtrap(x)=1000000*-(Uxz(Imax,x+1)-Uxz(Imax,x-1))/(2*dx);
end
for z=5:(gridsize-5)
        Fztrap(z)=-(Uxz(z+1,(gridsize/2))-Uxz(z-1,(gridsize/2)))/(2*dx);
end

Fgrav=Fg(1,:);
figure(7), clf(7), hold on
plot(Xtrap,Fxtrap,'-r','LineWidth',2)
%plot([-0.413,-0.413],[-10,10],'--k','linewidth',2) %68.98
%plot([0.413,0.413],[-10,10],'--k','linewidth',2) %68.98
%plot(Xtrap,Fgrav,'-k')
%plot(Xtrap,-Fgrav,'-k')
%plot(0.4255,-68.98,'-xb','linewidth',2)
xlabel('x [mm]')
%axis([-6,6,-10,10])
grid on
ylabel('F_{rad,lat} [\muN]')
hold off
%Fgrav(1)

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

