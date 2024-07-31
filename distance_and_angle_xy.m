function [d,theta] = distance_and_angle_xy(Xloc,Yloc,Zplane,xarray,yarray,zarray,Radius)
 d=sqrt((Zplane-zarray)^2+(Yloc-yarray)^2+(Xloc-xarray)^2);
 theta=acos((1/(Radius*d))*((Xloc-xarray)*(0-xarray)+(Yloc-yarray)*(0-yarray)+(Zplane-zarray)*(Radius-zarray)));
end