function [impar2d,ip_mean2d,ip_std2d,impar3d,ip_mean3d,ip_std3d,imparz,...
    ip_meanz,ip_stdz,deltazaxis]=impact_par(tracknr,radius,vertex,param_fit,Cf,N)

%   IMPACT_PAR calculated the impact parameters for every single track at
%   the inner side of the beam tube
%   First the the projected (x-y) impact parameters are calculated, then
%   the three dimensional impact parameter using first order approximation.
%   (tangent to track)
%
%   Input:  Radius          Radius of the beam tube
%           vertex          Cartesian coordinates of the vertices
%           param_fit		Fitted parameters at the inner side of the beam
%							tube
%           tracknr         array that contains the number of tracks per event
%                           that were used to calculate MC-pulls
%
%   Output:	impar2d			Array that contain the impact parameters for
%							all tracks projected to the x-y plane
%			ip_mean2d		Mean value of impar2d
%			ip_std2d		Standard deviation of impar2d
%           impar3d         Array that contains the 3-dimensional impact
%                           parameters
%           ip_mean3d       Mean value of impar3d
%           ip_std3d        Standard deviation of impar3d
%           imparz          Array that contains the z-component of the 3d 
%                           impact parameters
%           ip_meanz        Mean value of imparz
%           ip_stdz         Standart deviation of imparz
    
Phi_c	=  param_fit(:,1);
z_c		=  param_fit(:,2);
theta_c	=  param_fit(:,3);
beta_c	=  param_fit(:,4);
kappa_c	=  param_fit(:,5);

% no correction for 2*pi errors since phi_c goes only into 
% trigonometric functions
phi_c = beta_c+Phi_c;

% start parameters
x_c=radius*cos(Phi_c);
y_c=radius*sin(Phi_c);

x0=zeros(size(x_c));
y0=zeros(size(y_c));
z0=zeros(size(y_c));

% setting up an arry with the size of x_c with 
% the vertex coordinates for all tracks
index=1;
for eventnr=1:N.Event
    for k0=1:tracknr(eventnr)
		x0(index) = vertex(eventnr,1);
		y0(index) = vertex(eventnr,2);
		z0(index) = vertex(eventnr,3);
		index=index+1;
    end
end

hrH=1./kappa_c;
h=sign(hrH);

x_H=x_c-hrH.*sin(phi_c);
y_H=y_c+hrH.*cos(phi_c);

d = sqrt((x_H-x0).^2 + (y_H-y0).^2);
impar2d= hrH-h.*d;
%absimpar2d=abs(impar2d);
ip_mean2d = mean(impar2d);
ip_std2d  = std(impar2d);

% 3-dimensional impact parameter

Phi_p=atan2(y0-y_H,x0-x_H);
phi_p=Phi_p+h.*pi/2;

delta_phi = [phi_p-phi_c-4*pi,phi_p-phi_c-2*pi,phi_p-phi_c,...
             phi_p-phi_c+2*pi,phi_p-phi_c+2*2*pi];

dphi=zeros(size(phi_c));

for l=1:length(phi_c)
    k=find( delta_phi(l,:) >= -pi & delta_phi(l,:) <=pi );
    dphi(l)=delta_phi(l,k);
end


% l projected path length (=s*sin(theta))
% l= hrH.*dphi;
x_p = x0 + impar2d.*sin(phi_p); %+ l.*cos(phi_p) - (l.^2./(2*hrH)).*sin(phi_p);
y_p = y0 - impar2d.*cos(phi_p); %+ l.*sin(phi_p) + (l.^2./(2*hrH)).*cos(phi_p);
z_p = z_c + hrH.*cot(theta_c).*dphi;
% (x_p,y_p,z_p) is the point on the helix corresponding to the projected
% impact point (in general z_p~=z_0!)

% tangent in (x_p,y_p,z_p)
t_x=sin(theta_c).*cos(phi_p); 
t_y=sin(theta_c).*sin(phi_p);
t_z=cos(theta_c); 

s=sum([t_x t_y t_z].*[x_p-x0 y_p-y0 z_p-z0],2);
x_s=x_p-s.*t_x;
y_s=y_p-s.*t_y;
z_s=z_p-s.*t_z;

impar3d= sqrt((x_s-x0).^2+(y_s-y0).^2+(z_s-z0).^2);
%impar3d= sqrt((x_s-x0).^2+(y_s-y0).^2) *sign(-0.5+rand(1));

ip_mean3d=mean(impar3d,1);
ip_std3d=std(impar3d,[],1);
imparz= abs(z_p-z0);
ip_meanz=mean(imparz,1);
ip_stdz=std(imparz,[],1);


% 2nd interation - only for testing purpose in the forward region
% 

%%%%%%%%%% Estimation of delta(z) at the axis %%%%%%%%%%%%%
[N,M,O]=size(Cf);
dzbt=reshape(Cf(:,2,2),1,N);
dthetabt=reshape(Cf(:,3,3),1,N);
theta=reshape(param_fit(:,3),1,N);
deltazaxis=sqrt(dzbt.^2+radius^2*[1+1./(tan(theta).^2)].^2.*dthetabt.^2);
