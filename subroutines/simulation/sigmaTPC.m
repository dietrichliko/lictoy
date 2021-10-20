function sigsmear=sigmaTPC(paramf,zmin,zmax,k)

global convf sig TPC
global ITrack IEvent
global unit hwait

% Get neccesary data from input
z=paramf(2);
theta=paramf(3);
beta=paramf(4);
sigRPhi0=sig.RPhi0(k);
sigRPhi1=sig.RPhi1(k);
sigz0=sig.z0(k);
%sigz1=sig.z1(k);
CdiffRPhi=sig.CdiffRPhi(k);
Cdiffz=sig.Cdiffz(k);

% calculate missing parameters
if paramf(2)>0
    Ldrift=abs(z-zmax);
else
    Ldrift=abs(z-zmin);
end
Ldrift=Ldrift*1e-3*unit; % Ldrift in [m] for diffusion formula

if TPC.Number>1
    h=(max(TPC.Radius)-min(TPC.Radius))/TPC.Number; % padrow height
else
    h=6;
end
h=h/unit;               % h in [mm]

sigsmear(1)=sqrt( sigRPhi0^2 + sigRPhi1^2*sin(beta)^2 + CdiffRPhi^2 * (6/h) * sin(theta) * Ldrift);
sigsmear(2)=sqrt( sigz0^2 + Cdiffz^2 * Ldrift);