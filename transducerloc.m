function [xarray,yarray,zarray,phi,Nt] = transducerloc(NR,AR,Zi,Rot,M,Trap,deg1,deg2,deg3,deg4,deg5,deg6,deg7,deg8,deg9,deg10)

%Location and phase of the center transducer (origin)
xarray(1)=[0];
yarray(1)=[0];
zarray(1)=[0];
phi(1)=[0];
Nt=length(xarray);    %Counting transducers

%Location and phase of the 1st ring
for n=(Nt+1):(Nt+NR(1))
    xarray(n)=AR(1)*cos(2*pi*n/NR(1)+Rot(1)*pi/NR(1)+(0*35/180)*pi);%145/180)*pi
    yarray(n)=AR(1)*sin(2*pi*n/NR(1)+Rot(1)*pi/NR(1)+(0*35/180)*pi);%145/180)*pi
    zarray(n)=Zi(1);
    phi(n)=2*M*pi*n/NR(1)+Rot(1)*pi/NR(1);
    Nt=Nt+1;
end
%Location and phase of the 2nd ring
for n=(Nt+1):(Nt+NR(2))
    xarray(n)=AR(2)*cos(2*pi*n/NR(2)+Rot(2)*pi/NR(2));
    yarray(n)=AR(2)*sin(2*pi*n/NR(2)+Rot(2)*pi/NR(2));
    zarray(n)=Zi(2);
    phi(n)=2*M*pi*n/NR(2)+Rot(2)*pi/NR(2);
    Nt=Nt+1;
end
%Location and phase of the 3rd ring
for n=(Nt+1):(Nt+NR(3))
    xarray(n)=AR(3)*cos(2*pi*n/NR(3)+Rot(3)*pi/NR(3));
    yarray(n)=AR(3)*sin(2*pi*n/NR(3)+Rot(3)*pi/NR(3));
    zarray(n)=Zi(3);
    phi(n)=2*M*pi*n/NR(3)+Rot(3)*pi/NR(3);
    Nt=Nt+1;
end
%Location and phase of the 4th ring
for n=(Nt+1):(Nt+NR(4))
    xarray(n)=AR(4)*cos(2*pi*n/NR(4)+Rot(4)*pi/NR(4));
    yarray(n)=AR(4)*sin(2*pi*n/NR(4)+Rot(4)*pi/NR(4));
    zarray(n)=Zi(4);
    phi(n)=2*M*pi*n/NR(4)+Rot(4)*pi/NR(4);
    Nt=Nt+1;
end
switch Trap
    case 1
        phi=phi;
    case 2
        Phasedeg1=deg1;
        Phasedeg2=deg2;
        Phasedeg3=deg3;
        Phasedeg4=deg4;
        Phasedeg5=deg5;
        Phasedeg6=deg6;
        Phasedeg7=deg7;
        Phasedeg8=deg8;
        Phasedeg9=deg9;
        Phasedeg10=deg10;
        PhaseAction1=(2*pi*Phasedeg1)/360;
        PhaseAction2=(2*pi*Phasedeg2)/360;  
        PhaseAction3=(2*pi*Phasedeg3)/360;
        PhaseAction4=(2*pi*Phasedeg4)/360;
        PhaseAction5=(2*pi*Phasedeg5)/360;
        PhaseAction6=(2*pi*Phasedeg6)/360;
        PhaseAction7=(2*pi*Phasedeg7)/360;
        PhaseAction8=(2*pi*Phasedeg8)/360;
        PhaseAction9=(2*pi*Phasedeg9)/360;
        PhaseAction10=(2*pi*Phasedeg10)/360;
        phi(2)=PhaseAction1;
        phi(3)=PhaseAction2;
        phi(4)=PhaseAction3;
        phi(5)=PhaseAction4;
        phi(6)=PhaseAction5;
        phi(7)=PhaseAction6;
        phi(8)=PhaseAction7;
        phi(9)=PhaseAction8;
        phi(10)=PhaseAction9;
        phi(11)=PhaseAction10;
    case 3                      % Twin trap for Ultraino
        Nt=1;
        offset = pi/4;
        for i = 1:4 % # of ring
            
            if i == 1
                hr = NR(i)/2; % half ring value
            else
                if i == 2
                    hr = NR(i-1)+NR(i)/2; % half ring value
                else
                    if i == 3
                    hr = NR(i-2) + NR(i-1)+NR(i)/2; % half ring value
                    else
                    hr = NR(i-3) + NR(i-2)+ NR(i-1)+NR(i)/2; % half ring value
                    end
                end
            end
                
            for n=(Nt+1):(Nt+NR(i))
                if n < hr+2
                phi(n)=pi*0+Rot(i)*pi/NR(i);
                else
                phi(n)=pi*1+Rot(i)*pi/NR(i);
                end
                Nt=Nt+1;
            end
        end
      
end
end