function [qaramf,Cq]=vertconv(Bz,R,parami,Cp,convf,unit)

% function VERTCONV
% Called by vertexout
% Main program: LDT_main
%
% Input:    Bz      z value of magnetic field
%           R       Radius, where input parameters are defined
%           parami  Parameters (Phi,z,theta,beta,kappa) to be converted
%           Cp      Corresponding 5x5 covariance to be converted
%
% Output:   qaramf  Converted cartesian parameters (x,y,z,px,py,pz)
%           Cq      Converted 6x6 covariance matrix
%
% VERTCONV converts DELPHI like track parameters (Phi,z,theta,beta,kappa),
% defined at a cylindric reference surface of radius R, to cartesian
% parameters (x,y,z,px,py,pz). Also the corresponding 5x5 covariance matrix
% of the input parameters is converted to the 6x6 covariance matrix of rank
% 5, that corresponds to the converted cartesian parameters.

%global convf  % conversion factor [Gev/c T^(-1) m^(-1)]
%global unit

% readout of input parameters
Phi=parami(1);
z=parami(2);
theta=parami(3);
beta=parami(4);
kappa=parami(5);
phi=beta+Phi;
Q=-sign(kappa*Bz);      % charge
pt=-convf*Q*Bz/kappa;   % transverse momentum

% coordinates
q1=R*cos(Phi)/10*unit;  % [cm]
q2=R*sin(Phi)/10*unit;  % [cm]
q3=z/10*unit;           % [cm]

% momentum
fac=convf*abs(Bz)/abs(kappa);
px=cos(phi)*fac;
py=sin(phi)*fac;
pz=cot(theta)*fac;

% derivative matrix elements
dxPhi=-q2;
dyPhi=q1;
dzz=unit/10;				 
% and for the momenta
dpxPhi=-py;
dpyPhi=px;
dpztheta=-1/(sin(theta))^2*fac;
dpxbeta=-py;
dpybeta=px;	
dpxkappa=-px/kappa;
dpykappa=-py/kappa;
dpzkappa=-pz/kappa;

% derivative matrix
D=zeros(6,5);

D(1,1)=dxPhi;
D(2,1)=dyPhi;
D(4,1)=dpxPhi;
D(5,1)=dpyPhi;
D(3,2)=dzz;
D(6,3)=dpztheta;
D(4,4)=dpxbeta;
D(5,4)=dpybeta;
D(4,5)=dpxkappa;
D(5,5)=dpykappa;
D(6,5)=dpzkappa;

% converted covariance matrix
Cq=D*Cp*D';

% converted parameters
qaramf=[q1,q2,q3,px,py,pz];