function [paramf,dermat,alrphi,ierr]=plan5mod(zi,parami,zf,rmin,rmax,iopt)
%                                                                      *
%                                                                      *
%    aim :                                                             *
%    -----                                                             *
%    extrapolate a helix defined by the initial parameters parami      *
%    up to a given z-plane, and compute the derivatives of the         *
%    final parameters w.r.t. the initial ones.                         *
%    Note: no backward tracking!                                       *
%                                                                      *
%    the computation uses double precision on intermediate variables   *
%    if the variation of phi angle is less than dphimn (.0001 in this  *
%    version) the computation is done at first order in 1/r in order   *
%    to avoid rounding errors in the derivatives                       *
%                                                                      *
%    input  :  zi            : initial z                               *
%    input  :  parami(1:5)   : initial parameters                      *
%                              (x,y,theta,phi,1/r signed)              *
%              zf            : z of the final plane                    *
%              rmin          : lower limit of r on the plane           *
%              rmax          : upper limit of r on the plane           *
%              iopt          : 0 if derivatives not requested          *
%                              1 if derivatives requested              *
%                                                                      *
%    output :  ierr          : 0 if ok                                 *
%                              1 if no intersection found              *
%                              3 if intersection outside of limits     *
%              paramf(1:5)   : final parameters                        *
%              dermat(5,5)   : deriv. of final w.r.t. initial param.   *
%                              der(1) = d(x)/d(theta)                  *
%                              der(2) = d(x)/d(phi)                    *
%                              der(3) = d(x)/d(1/r)                    *
%                              der(4) = d(y)/d(theta)                  *
%                              der(5) = d(y)/d(phi)                    *
%                              der(6) = d(y)/d(1/r)                    *
%                              der(7) = d(phi)/d(theta)                *
%                              der(8) = d(phi)/d(1/r)                  *
%              alrphi        : length (in r-phi projection) from start *
%                              to extrapolation, with a sign (positive *
%                              if the extrapolation is towards the     *
%                              direction defined by theta,phi)         *
%                                                                      *
%    author  :  p. billoir                                             *
%                                                                      *
%    first version : 26-01-88                                          *
%                                                                      *
%
twopi=pi*2;
dphimn=1.0e-04;
dermat=eye(5);
%
ierr=0;
%
dz=zf-zi;
tanth=tan(parami(3));
rdphi=dz*tanth;
%
alrphi=rdphi;
%
%   x,y coordinates of intersection
rtrk=1./parami(5);
cosf0=cos(parami(4));
sinf0=sin(parami(4));
xc=parami(1)-rtrk*sinf0;
yc=parami(2)+rtrk*cosf0;
dphi=parami(5)*rdphi;
phi1=mod(parami(4)+dphi,twopi);
cosf1=cos(phi1);
sinf1=sin(phi1);
x1=xc+rtrk*sinf1;
y1=yc-rtrk*cosf1;
r1=sqrt(x1^2+y1^2);
%
%   intersection outside of limits in r
if (r1<rmin|r1>rmax)
    ierr=3;
end
%   parameters at the intersection
paramf(1)=x1;
paramf(2)=y1;
paramf(3)=parami(3);
paramf(4)=phi1;
if (paramf(4)<0.)
    paramf(4)=paramf(4)+twopi;
end
paramf(5)=parami(5);
%
%   computation of derivatives -----------------------------------
%
if (iopt==1)
    ct2inv=1.+tanth^2;
    der(1)=ct2inv*dz*cosf1;
    der(4)=ct2inv*dz*sinf1;
    der(7)=dz*parami(5)*ct2inv;
    der(8)=rdphi;
    %
    %   "exact" formulae if |dphi| > dphimn
    if (abs(dphi)>=dphimn)
        dcosf=cosf1-cosf0;
        dsinf=sinf1-sinf0;
        der(2)=rtrk*dcosf;
        der(3)=rtrk^2*(dphi*cosf1-dsinf);
        der(5)=rtrk*dsinf;
        der(6)=rtrk^2*(dphi*sinf1+dcosf);
        %
        %   first order in 1/r if |dphi| < dphimn
    else
        der(2)=-rdphi*sinf0;
        der(3)=.5*rdphi*der(2);
        der(5)=rdphi*cosf0;
        der(6)=.5*rdphi*der(5);
    end
    dermat(1,3)=der(1);
    dermat(1,4)=der(2);
    dermat(1,5)=der(3);
    dermat(2,3)=der(4);
    dermat(2,4)=der(5);
    dermat(2,5)=der(6);
    dermat(4,3)=der(7);
    dermat(4,5)=der(8);
    %der(1) = d(x)/d(theta)                  *
    %der(2) = d(x)/d(phi)                    *
    %der(3) = d(x)/d(1/r)                    *
    %der(4) = d(y)/d(theta)                  *
    %der(5) = d(y)/d(phi)                    *
    %der(6) = d(y)/d(1/r)                    *
    %der(7) = d(phi)/d(theta)                *
    %der(8) = d(phi)/d(1/r)                  *
end