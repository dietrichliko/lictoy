function [parf,z,dermat]=cyl2plan(pari,R,ider)
%
% input:    pari(1:5)   :   Phi, z, theta, beta, kappa (=1/R)
%           R           :   Radius of zylinder layer
%           ider        :   1: derivative matrix requested, 0: not
% output:   parf        :   x, y, theta, psi, kappa
%           z           :   z value, where the new forward parameters are
%                           defined
%           dermat      :   derivative matrix of conversion
%
% CYL2PLAN converts a given parameter set from a "cylinder" representation
% to a forward ("plane") representation  
% (Phi,z,theta,beta,kappa) at R ==> (x,y,theta,phi,kappa) at z
%
% author: M.Regler [1]
% [1] R. Fruehwirth, P.K. Lichtenwagner, M. Regler and D. Stampfer: 'the
% delphi forward track fit; track fitting with outlier rejection', NIM in
% Physics A334 (1993) 528-536


% This file is the same as cylpane.m. It is used to change the derivative 
% matrix independently at the barrel to plane and plane to TPC region
dermat=eye(5);
z=pari(2);

tanth=tan(pari(3));
Phi=pari(1);
%Phi=RPhi/R;
psi=pari(4)+Phi;
sinpsi=sin(psi);
cospsi=cos(psi);
sinPhi=sin(Phi);
cosPhi=cos(Phi);
x=R*cosPhi;
y=R*sinPhi;

parf(1)=x;
parf(2)=y;
parf(3)=pari(3);
parf(4)=psi;
parf(5)=pari(5);

%----Derivative matrix----

dermat(1,1)=-sinPhi*R;
dermat(1,2)=-tanth*cospsi;
dermat(2,1)=cosPhi*R;
dermat(2,2)=-tanth*sinpsi;
dermat(4,1)=1;
dermat(4,2)=-tanth*pari(5);

ierr=0;
