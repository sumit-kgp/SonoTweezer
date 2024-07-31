function [d,theta] = distance_and_angle_yz(Xaxis,Yloc,Zloc,xarray,yarray,zarray,Focus)
 d=sqrt((Zloc-zarray)^2+(Yloc-yarray)^2+(Xaxis-xarray)^2);
 theta=acos((1/(Focus*d))*((Xaxis-xarray)*(0-xarray)+(Yloc-yarray)*(0-yarray)+(Zloc-zarray)*(Focus-zarray)));
end
% 
% function [d,theta] = distance_and_angle_yz(Xloc,Yloc,Zplane,xarray,yarray,zarray,Radius)
%  d=sqrt((Zplane-zarray)^2+(Yloc-yarray)^2+(Xloc-xarray)^2);
%  theta=acos((1/(Radius*d))*((Xloc-xarray)*(0-xarray)+(Yloc-yarray)*(0-yarray)+(Zplane-zarray)*(Radius-zarray)));
% end