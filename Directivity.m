clear all, clc, close all
%%
%Initiation
%
%Constants
c0air=343;            %Speed of sound in air        [m/s]
c0water=1510
fa=40000;              %Transducer frequency         [Hz]
fw=1000000;
a=0.0045;             %Piston radius                [m]
a2=0.0065;

%%
%Calculation of directivity function for polarplot
%
lambda_air=c0air/fa;              %Wavelength in air              [m]
lambda_water=c0water/fw;
k_air=(2*pi)/lambda_air;         %Wave number for 40kHz in air   [1/m]
k_water=(2*pi)/lambda_water;                    %Wave number for 80kHz in air   [1/m]
thetad = linspace(-pi/2,pi/2,1000);  %Define grid for angle          [rad]

for t=1:length(thetad)
[DirectFunc_1(t)]=directfunc(k_air,a,thetad(t));
[DirectFunc_2(t)]=directfunc(k_water,a,thetad(t));
end
figure(1), clf(1), hold on
subplot(1,2,1)
polarplot(thetad,DirectFunc_1,'-r','LineWidth',2)
title('Air 40kHz')
subplot(1,2,2)
polarplot(thetad,DirectFunc_2,'-b','LineWidth',2)
title('Water 1.0MHz')
hold off


