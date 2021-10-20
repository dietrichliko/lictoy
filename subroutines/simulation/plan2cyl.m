function [paramcyl,R,Dpc]=plan2cyl(parami,z0,ider)

% input:    parami      :   x,y,theta, phi, kappa (=1/R signed)
%          
% output:   parmcyl     :   Phi, z, theta, beta, kappa
%           R           :   2 dimensional radius 
%           D           :   derivative matrix of conversion
%
% PLAN2CYL converts a given parameter set from a "plane" representation
% to a cylindrical representation  
% (x,y,theta,phi,kappa) at z0 ==> (Phi,z,theta,beta,kappa) at R

Dpc=zeros(5);

% to make the equations easier to read
x       =   parami(1);
y       =   parami(2);
z       =   z0;
phi     =   parami(4);
R       =   sqrt(x^2+y^2);
theta   =   parami(3);
%theta = atan(R/z0);

% new coordinates
Phi=atan2(y,x);
beta=phi-Phi;            

% forcing beta to the interval [-pi/2,pi/2]
if beta>pi/2
    beta=beta-2*pi;
end
if beta<-pi/2
    beta=beta+2*pi;
end

% ---------------------------- DERIVATIVE MATRIX --------------------------
if ider
    % first compute derivative matrix Dcp of conversion cylinder->plane
    % and invert it afterwards
    Dcp=eye(5);
    Dcp(1,1)=-sin(Phi)*R;
    Dcp(1,2)=-tan(theta)*cos(phi);
    Dcp(2,1)=cos(Phi)*R;
    Dcp(2,2)=-tan(theta)*sin(phi);
    Dcp(4,1)=1;
    Dcp(4,2)=-tan(theta)*parami(5);
    
    Dpc=inv(Dcp);

    %{
    D(1,1)=-y/R^2;
    D(1,2)=x/R^2;
    D(2,1)=x/R*cot(parami(3))*cos(beta);
    D(2,2)=y/R*cot(parami(3))*cos(beta);
    D(2,3)=-R/(sin(parami(3)))^2*cos(beta);
    D(2,4)=-R*cot(parami(3))*sin(beta);
    D(3,3)=1;
    D(4,1)=y/R^2;
    D(4,2)=-x/R^2;
    D(4,4)=1;
    D(5,5)=1;
    %}
end

paramcyl=[Phi,z,parami(3),beta,parami(5)];
