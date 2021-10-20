function [paramr,Hk,Vk,Ak,varMS]=reftrack(pstartMS,param_start,reftype,...
                                                  refindex,ref)

% REFTRACK 
%
% function [paramr,Hk,Vk,Ak,varMS,measnr1]=REFTRACK(pstartMS, param_start,
%                                                    reftype,refindex,ref)
% Called by LTD_main
% Main program: LTD_main
%
% Input:    pstartMS    Start parameters at outer side of beamtube (incl.MS)
%                       (Phi,z,theta,beta,kappa)
%           param_start Start parameters at inner side of beamtube (excl.MS)
%                       (Phi,z,theta,beta,kappa)
%           reftype     Array determining the type of reference surface
%                       according to the simulation:
%                           1: cylinder layer
%                           0: plane forward layer
%           refindex    holds the index of the corresponding cylinder layer
%                       or forward layer for each reference surface, e.g if
%                       the 4th layer traversed is the 2nd forward layer,
%                        -> refindex(4)=2.
%           ref         holds the corresponding geometrical value of the
%                       reference surfaces. In case of a barrel layer it
%                       holds its radius, in case of a forward layer it
%                       holds its z position
%
% Output:   paramr      Array of reference parameters, defined at each
%                       reference surface, in the corresponding
%                       representation. (Phi,z,theta,beta,kappa) for a
%                       cylinder layer and (x,y,theta,phi,kappa) for a
%                       forward plane layer
%           Hk          Array of H matrices for every layer
%           Vk          Array of measurement error matrices
%           Ak          Array of derivative matrices for every layer
%           varMS       Variances of multiple scattering
%
% REFTRACK computes the reference track (the expansion points for the 
% Kalman filter), the H matrices, the error matrices, the derivative
% matrices and the variances of multiple scattering for every layer hit
% according to the simulation.
% The reference track is the completely undisturbed extrapolation of the
% start parameters through the whole detector, without kinks due to 
% multiple scattering.

global Bz convf sinbmx unit RadiusIR ITrack Mass numhitnr1 hwait Flags 
global fidlog whandle mhandle N radius zpos delta Xlen detnr1 eff distr 
global sig d delta unit name

warning off;

% prepare data of forward/rear region, respectively
if param_start(3)>pi/2 % rear direction
    effu=eff.ru;effv=eff.rv;deltau=delta.ru;deltav=delta.rv;fXlen=Xlen.r;
    fdistr=distr.r;FLayer=N.RLayer;sigu=sig.ru;sigv=sig.rv;du=d.ru;dv=d.rv;
else
    effu=eff.fu;effv=eff.fv;deltau=delta.fu;deltav=delta.fv;fXlen=Xlen.f;
    fdistr=distr.f;FLayer=N.FLayer;sigu=sig.fu;sigv=sig.fv;du=d.fu;dv=d.fv;
end
% prepare data of barrel region
Radius=radius.b;bdistr=distr.b;bXlen=Xlen.b;BLayer=N.BLayer;

% specifications for the helix propagation
iopt=1;     % derivative matrix requested
idir=1;     % forward propagation
rmax=inf;   % real dimensions do not matter, as we already know from the 
rmin=0;     % simulation, which layers have been hit
%zmax=zpos.bmax(refindex(reftype));
%zmin=zpos.bmin(refindex(reftype));

% special treatment of the beamtube
%pstartMS=param_start;
parami=pstartMS;         % initial parameters are the params. at the BT
%paramr(1,:)=param_start;    % first ref. params. are the params. at the BT
paramr(1,:)=pstartMS;    % first ref. params. are the params. at the BT
Hk(1,:,:)=zeros(2,5)*NaN;   % no measurement at BT -> no H matrix
Vk(1,:,:)=diag([NaN,NaN]);  % no measurement at BT -> no error matrix
Ak(1,:,:)=eye(5);           % first derivative matrix is the unit matrix, 
                            % it just fills the space, since there is no
                            % further propagation step when having reached
                            % the beam tube at filtering inwards in the
                            % Kalman filter


pr=abs(convf*Bz/paramr(1,5)/sin(paramr(1,3)));  % Momentum of reference track
measNr1=length(ref);                            % outermost measurement

%-- prepare an array to determine what to do: conversions, propagations --%
reftype1=[reftype,reftype(end)];    % shift to the left
reftype1(1)=[];
change=reftype1-reftype;            % ~=0 when conversion needed
prop=reftype1+reftype;              % ~=1 when normal propagation
warning off
L=logical(change);
%warning on
prop(L)=change(L);                  % writes kind of conversion
prop(end)=[];                       % n layers -> n-1 propagations!
% prop determines the action needed to do:
%   -1: conversion barrel-forward and forward propagation
%    0: no conversion needed, only forward propagation
%    1: conversion forward-barrel and barrel propagation
%    2: no conversion needed, only barrel propagation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:length(ref)     % loop over all layers
    % find out what to do using prop
   switch prop(k-1)
       case -1
           % last layer was barrel layer, next one is forward layer
           % conversion barrel-forward and forward propagation
           [parami,zi,Aconv]=cyl2plan(parami,ref(k-1),iopt);
           bprop=0;
       case 0
           % both last and next layer are forward layers
           % no conversion needed, only forward propagation
           Aconv=eye(5);
           zi=ref(k-1);
           bprop=0;
       case 1
           % last layer was forward layer, next one is barrel layer
           % conversion forward-barrel and barrel propagation
           [parami,Ri,Aconv]=plan2cyl(parami,ref(k-1),iopt);
           bprop=1;
       case 2
           % both last and next layers are barrel layers
           % no conversion needed, only barrel propagation
           Aconv=eye(5);
           Ri=ref(k-1);
           bprop=1;
   end % switch prop(k-1)
   
   % now conversions are done, propagation follows
   %----------------------------------------------------------------------%
   if bprop
       % barrel propagation, H matrix and error matrix
       
       % parameters are prepared, normal cylinder propagation
       [paramf,Aprop,path,ierr]=prop5mod(Ri,parami,...
                    idir,ref(k),-inf,inf,sinbmx,iopt);
        
        % prepare data for H matrix
        bnow=refindex(k);
        R=Radius(bnow);
        alpha=delta.b(bnow);
        Hk(k,:,:)=[R            0          0 0 0   % H matrix for this
                   R*cos(alpha) sin(alpha) 0 0 0]; % layer including alpha
               
        switch bdistr(bnow)
            case 0 % simulate gaussian distributed errors
                sigsmear=[sig.RPhi0(bnow),sig.z0(bnow)];
                % random measurement errors according to computed sigmas
            case 1    % simulate uniformly distributed errors
                % random measurement errors
                sigsmear(1)=d.RPhi(bnow)/sqrt(12); % computes equivalent sigmas
                sigsmear(2)=d.z(bnow)/sqrt(12);    % from strip distances
            case 2 % simulate special TPC errors
                sigsmear=sigmaTPC(paramf,zpos.bmin(bnow),zpos.bmax(bnow),bnow);
            otherwise
                sigsmear=[NaN,NaN];
        end % if uniform(k)==0 - end
        
        % end of barrel part
   %----------------------------------------------------------------------%
   else % if bprop
       % forward propagation, H matrix and error matrix
       
       % parameters are prepared, normal forward plane propagation
       [paramf,Aprop,path,ierr]=plan5mod(zi,parami,...
           ref(k),rmin,rmax,iopt);
       fnow=refindex(k);
       
       % Derive H-matrix
       % Step 1: derive polar angle Phi
       Phi=atan(paramf(2)/paramf(1));
       if paramf(1)<0       Phi=pi+Phi;
       elseif paramf(2)<0   Phi=2*pi+Phi; end

       % Step 2: add angles defining the directions
       Psi1=Phi+deltau(fnow); Psi2=Phi+deltav(fnow);
       % make sure that both are between 0 and 2*pi
       while Psi1>2*pi Psi1=Psi1-2*pi; end
       while Psi1<0    Psi1=Psi1+2*pi; end
       while Psi2>2*pi Psi2=Psi2-2*pi; end
       while Psi2<0    Psi2=Psi2+2*pi; end

       % H matrix of current layer, dependent on intersection point
       Hk(k,:,:)=[cos(Psi1) sin(Psi1) 0 0 0
                  cos(Psi2) sin(Psi2) 0 0 0];

       if fdistr(fnow)==0    % simulate gaussian distributed errors
           sigsmear(1)=sigu(fnow);
           sigsmear(2)=sigv(fnow);
       elseif fdistr(fnow)==1    % simulate uniformly distributed errors
           sigsmear(1)=du(fnow)/sqrt(12);    % computes equivalent sigmas
           sigsmear(2)=dv(fnow)/sqrt(12);    % from strip distances
       else
           sigsmear=[NaN,NaN];
       end % if fdistr(refindex(k))==0

       % end of forward part
              
   end  % if bprop
   %----------------------------------------------------------------------%
   
   % last steps inside loop over the layers
   if ierr
       % if at any time the propagation fails, the rest of
       break;     % the reference track is obsolete and will be ignored
   end
   
   Ak(k,:,:)=Aprop*Aconv; % stores the deriv. matrix
   Vk(k,:,:)=diag(sigsmear.^2);       % stores the error matrix
   paramr(k,:)=paramf;                % stores the reference parameters
   parami=paramf;                     % new parameters for next loop

end % for k=2:length(ref) % loop over all layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare variances of multiple scattering
if Flags.MulSca~=0
    b=reftype(1:size(paramr,1)); % first radiation length of barrel layers
    if sum(b)~=0
      X(b)=bXlen(refindex(b))./(sin(paramr(b,3)).*cos(paramr(b,4)))';
      % Effective thickness of scatterer incl. beta (all layers)
    end
    f=~b;   % now radiation length of forward layers
    if sum(f)~=0
      X(f)=fXlen(refindex(f))./abs(cos(paramr(f,3)))';
    end
    
    sigMS=0.0136*sqrt((Mass^2+pr^2)/pr^4)*sqrt(X).*(1+0.038*log(X));
          % s.d. of projected multiple scattering angle
    varMS=sigMS.^2;
else
    varMS=NaN;
end % if Flags.MulSca~=0