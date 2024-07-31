function [d,theta] = distance_and_angle(Xloc,Yaxis,Zloc,xarray,yarray,zarray,Focus)
 d=sqrt((Zloc-zarray)^2+(Yaxis-yarray)^2+(Xloc-xarray)^2);
 theta=acos((1/(Focus*d))*((Xloc-xarray)*(0-xarray)+(Yaxis-yarray)*(0-yarray)+(Zloc-zarray)*(Focus-zarray)));
end