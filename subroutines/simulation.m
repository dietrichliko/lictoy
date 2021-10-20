function [pstartMS,param_true,coord1,coord2,meas,hitpattern,...
          reftype,refindex,ref,Xstore,varMS_store,V,H]=simulation(param_start,reg)
      
% SIMULATION
%
% function [paramstart_MulSca,param_true,coord1,coord2,meas,hitpattern,...
%          reftype,refindex,ref]=SIMULATION(param_start,reg)
%
% Called by LDT_main
% Main program: LDT_main
%
% Input:    param_start Start parameters at the inner side of the beam tube
%                       (Phi,z,theta,beta,kappa)
%           reg         Flag determining the region the track will traverse
%                           ff: forward
%                           if: intermediate forward
%                           bb: barrel
%                           ir: intermediate rear
%                           rr: rear
%
% Output:   pstartMS            'Real' start parameters at the outer side of
%                               the beamtube
%           param_true          Parameters of true intersection points at 
%                               each layer, including multiple scattering,
%                               representation depending on the type of
%                               surface:    (Phi,z,theta,beta,kappa) for a
%                               cylinder layer and (x,y,theta,phi,kappa) 
%                               for a forward plane layer
%           coord1              Coord. 1 measurements (RPhi or u) for every
%                               layer, including multiple scattering and 
%                               measurement errors
%           coord2              Coord. 2 measurements (z or v) for every 
%                               layer, including multiple scattering and 
%                               measurement errors
%           meas                Indicates which layer measured which
%                               coordinate: 0: no measurement at all
%                                           1: coord1 measured
%                                           2: coord2 measured
%                                           3: both measured
%           hitpattern          holds the names of the layers in the order
%                               they have been traversed
%           reftype             Array determining the type of reference 
%                               surface according to the simulation:
%                                   1: cylinder layer
%                                   0: plane forward layer
%           refindex            holds the index of the corresponding 
%                               cylinder layer or forward layer for each 
%                               reference surface, e.g if the 4th layer 
%                               traversed is the 2nd forward layer, ->
%                                   refindex(4)=2.
%           ref                 holds the corresponding geometrical value 
%                               of the reference surfaces. In case of a 
%                               barrel layer it holds its radius, in case 
%                               of a forward layer it holds its z position
%
% SIMULATION simulates one track going through an arbitrary part of
% the chosen detector setup. It starts with randomly generated start
% parameters put in via param_start. These are the parameters at the inner 
% side of the beam tube. 
% (1) Multiple scattering is added according to the type of surface
% (2) The parameters are converted to both barrel and forward
%       representation
% (3) Both sets of parameters are propagated along a helix to the next
%       barrel or forward layer, respectively.
% (4) The propagation with the lower pathlength defines the next hit
%
% These steps are repeated until the track leaves the detector.
% In every step the measurements are computed from the simulated state
% vector via the Hk matrix and a random measurement error is added.
          
global Bz convf sinbmx unit RadiusIR ITrack IEvent Mass numhitnr1 hwait
global Flags fidlog whandle mhandle N radius zpos delta Xlen detnr1 eff
global delta unit name distr sig d

if param_start(3)>pi/2 % rear direction
    z=zpos.r;innerrad=radius.rin;outerrad=radius.rout;effu=eff.ru;
    effv=eff.rv;deltau=delta.ru;deltav=delta.rv;fXlen=Xlen.r;fdistr=distr.r;
    FLayer=N.RLayer;sigu=sig.ru;sigv=sig.rv;du=d.ru;dv=d.rv;
    fname=name.r;
else
    z=zpos.f;innerrad=radius.fin;outerrad=radius.fout;effu=eff.fu;
    effv=eff.fv;deltau=delta.fu;deltav=delta.fv;fXlen=Xlen.f;fdistr=distr.f;
    FLayer=N.FLayer;sigu=sig.fu;sigv=sig.fv;du=d.fu;dv=d.fv;
    fname=name.f;
end
Radius=radius.b;bname=name.b;bdistr=distr.b;bXlen=Xlen.b;BLayer=N.BLayer;

% specifications for the helix propagation
iopt=0;  % 0: derivative matrix NOT requested
idir=1;  % forward propagation

Xstore=NaN;
varMS_store=NaN;

% initial parameters are the parameters at the BT
parami=param_start;  
pstartMS=param_start;
param_true(1,:)=param_start;
hitpattern(1)=bname(1);
meas(1)=0;
coord1(1)=NaN;
coord2(1)=NaN;
ref(1)=Radius(1);   % stores the radius and the z value, respectively, of
                    % the corresponding reference surfaces
refindex(1)=1;      % stores the number of the corresponding ref. surfaces
reftype(1)=true;       % stores the type of the corresponding ref. surfaces
                    % 1 = cylinder surface, 0 = plane surface
bnow=1;             % present cylinder layer
fnow=1;             % present forward layer
Rnow=Radius(1);     % current radius
znow=param_start(2);% current z position
bprop=0;
fprop=0;
maxzmax=max(zpos.bmax);
minzmin=min(zpos.bmin);

for k=2:(FLayer+BLayer)     % loop over all layers, terminated by break statement
    if Flags.MulSca    
        pT=convf*Bz/parami(5); % Transv. momentum computed from curvature
        p=pT/sin(parami(3));   % absolute Momentum
        
        if reftype(k-1) X=bXlen(bnow)/(sin(parami(3))*cos(parami(4)));
        else            X=fXlen(fnow)/abs(cos(parami(3)));  end
        sigms=0.0136*sqrt((Mass^2+p^2)/p^4)*sqrt(X)*(1+0.038*log(X));
        % s.d. of projected multiple scattering angle
        Xstore(k-1)=X;  varMS_store(k-1)=sigms^2;
              
        ran=(randn(1,2).*[1 1/sin(parami(3))])*sigms; %changes in direction
        if parami(3)>pi/2 ran(1)=-ran(1); end
        
        
        parami(3)=parami(3)+ran(1); % Kick in theta
        parami(4)=parami(4)+ran(2); % Kick in phi (beta)
        pT=p*sin(parami(3));        % recompute pt with new theta
        parami(5)=convf*Bz/pT;      % Curvature from transverse momentum
        
        if k==2 pstartMS=parami; end
        
    end % if Flags.MulSca
    
    if ~reftype(k-1)
        % parami are in forward representation
        % convert to barrel representation for cyl. propagation
        [bparami,Ri,Dummy]=plan2cyl(parami,ref(k-1),iopt);
        fparami=parami;
        zi=z(fnow);
    else
        % parami are in barrel representation
        % convert to forward representation for plane propagation
        [fparami,zi,Dummy]=cyl2plan(parami,ref(k-1),iopt);
        bparami=parami;
        Ri=Radius(bnow);
    end

    % find out cylinder layers left to propagate to
    % only consider layers with R>Rnow
    bnext=min(find(Radius>Rnow));
    if ( ~isempty(bnext) & reg~=1 )
        % propagation to the next hit barrel layer
        %--------------------Layer search------------------------
        for bprop=bnext:BLayer
            zmin=zpos.bmin(bprop);
            zmax=zpos.bmax(bprop);
            Rf=Radius(bprop);
            [bparamf,Dummy,bpath,berr]=prop5mod(Ri,bparami,...
                idir,Rf,zmin,zmax,sinbmx,iopt);
            if berr~=3
                break;
            end
            
            if bparamf(2)>maxzmax | bparamf(2)<minzmin
                bpath=inf; break;
            end
            
            if bprop==BLayer bpath=inf; end
            % set pathlength to infinity if loop ends without having hit
            % any layer
        end % for bprop=bnext:BLayer - end
        %------------------------------------------------------------
        if berr~=0
            if ( (berr==2) & (bparamf(2)<=zmax & bparamf(2)>=zmin) )
            %if ( (berr==2 | berr==1) )
            %{
                str=get(whandle,'String');
                str{2}=['Partial Track ',num2str(ITrack),' Event ',num2str(IEvent)];
                str{3}=['sin(beta_max) exceeded at R=',num2str(Rf,'% 0.2f'),' Layer: ',bname{bprop}];
                str{4}=(['Reconstruction starts at ',bname{bprop-1}]);
                set(whandle,'String',str);
                set(whandle,'foregroundcolor','red')
                drawnow
                fprintf(fidlog,'%s\n',str{2});
                fprintf(fidlog,'%s\n',['   ',str{3}]);
                fprintf(fidlog,'%s\n',['   ',str{4}]);
            %}
            end
            bpath=inf;
            %break;
        end
    else
        % end of barrel region reached
        % or track in forward/rear region only
        bprop=BLayer;
        berr=1;
        bpath=inf;
    end

    % find out forward layers left to propagate to
    % only consider layers with |z|>|znow|
    fnext=min(find(abs(z)>abs(znow)));
    if ( ~isempty(fnext) & reg~=3 )
        % propagation to the next hit forward layer
        for fprop=fnext:FLayer
            rmin=innerrad(fprop);
            rmax=outerrad(fprop);
            zf=z(fprop);
            [fparamf,Dummy,fpath,ferr]=plan5mod(zi,fparami,...
                zf,rmin,rmax,iopt);
            if ferr==0 break; end % jump out of loop if a layer was hit
            if fprop==FLayer fpath=inf; end
            % set pathlength to infinity if loop ends without having hit
            % any layer
        end
    else
        % end of forward region reached
        % or track in barrel region only
        fprop=FLayer;
        ferr=1;
        fpath=inf;
    end

    %if ( bprop==BLayer & berr~=0 & fprop==FLayer & ferr~=0)
    %if ( bprop==BLayer & fprop==FLayer)
    if ( bpath==inf & fpath==inf )
        break;
    end

    % Now we have bparamf and fparamf, defined at the next hit barrel
    % and forward layer, respectively.
    % The propagation with the lower pathlength (bpath or fpath) is
    % the correct one to continue the simulation

    if bpath<fpath
        % next layer is a cylinder layer
        reftype(k)=true;
        bnow=bprop; % now we sit on the next barrel layer
        ref(k)=Radius(bnow);
        refindex(k)=bnow;
        hitpattern(k)=bname(bnow);
        param_true(k,:)=bparamf;
        paramf=bparamf;
        Rnow=Radius(bnow);
        znow=bparamf(2);

        % consider inefficiencies
        meas(k)=0;
        meas1=rand(1)<=eff.RPhi(bnow);
        if eff.z(bnow)==-1 meas2=meas1;
        else meas2=rand(1)<eff.z(bnow);
        end
        if meas1&meas2   meas(k)=3;
        else
            if meas2     meas(k)=2;
            elseif meas1 meas(k)=1;
            end
        end

        switch meas(k)
            case 0
                H0=zeros(2)*NaN;    % no measurement
            case 1
                H0=[1,0;0,0];   % only RPhi measured
            case 2
                H0=[0,0;0,1];   % only z measured
            otherwise
                H0=eye(2);      % both measured
        end
        % The matrix H0 eliminates measurements of inefficient layers or layers
        % that were not hit

        alpha=delta.b(bnow);
        Hk=[Radius(bnow)            0          0 0 0
            Radius(bnow)*cos(alpha) sin(alpha) 0 0 0];
        % H-Matrix incl. stereo angle alpha

        paramh=Hk*paramf;            % transformation vector space -> meas. sp.
        RPhi=paramh(1);              % local variables
        zb=paramh(2);
        meas_true=[RPhi;zb];

        switch bdistr(bnow)
            case 0 % simulate gaussian distributed errors
                
                sigsmear=[sig.RPhi0(bnow),sig.z0(bnow)];
                % random measurement errors according to computed sigmas
                smear=randn(1,2).*sigsmear;
                meas_store=meas_true+H0*smear'; % errors are added
                % H0 eliminates inefficient measurements by neglecting the
                % measurement error. So the 'measurement' of an inefficient
                % layer is the true intersection point. However, in the
                % reconstruction the measurements marked as failed do
                % not contribute!
                coord1(k)=meas_store(1);  % stores RPhi and z
                coord2(k)=meas_store(2);
                
            case 1    % simulate uniformly distributed errors

                % random measurement errors
                smear(1)=-d.RPhi(bnow)/2+rand(1)*d.RPhi(bnow);
                smear(2)=-d.z(bnow)/2+rand(1)*d.z(bnow);
                meas_store=meas_true+H0*smear'; % errors are added
                % concerning H0 see comment above
                coord1(k)=meas_store(1);  % stores RPhi and z
                coord2(k)=meas_store(2);
                sigsmear(1)=d.RPhi(bnow)/sqrt(12); % computes equivalent sigmas
                sigsmear(2)=d.z(bnow)/sqrt(12);    % from strip distances
                
            case 2 % simulate special TPC errors
                
                sigsmear=sigmaTPC(paramf,zmin,zmax,bnow);
                smear=randn(1,2).*sigsmear;
                meas_store=meas_true+H0*smear'; % errors are added
                coord1(k)=meas_store(1);  % stores RPhi and z
                coord2(k)=meas_store(2);
                
            otherwise
                
                sigsmear=[NaN,NaN];
                
        end % switch bdistr(bnow)
    else
        % next layer is a forward layer
        reftype(k)=false;
        fnow=fprop; % now we sit on the next forward layer
        ref(k)=z(fnow);
        refindex(k)=fnow;
        hitpattern(k)=fname(fnow);
        param_true(k,:)=fparamf;
        paramf=fparamf;
        Rnow=sqrt(fparamf(1)^2+fparamf(2)^2);
        znow=z(fnow);

        meas(k)=0;
        meas1=rand(1)<=effu(fnow);
        if effv(fnow)==-1   meas2=meas1;
        else                meas2=rand(1)<effv(fnow);
        end
        if meas1&meas2      meas(k)=3;
        else
            if meas2        meas(k)=2;
            elseif meas1    meas(k)=1;
            end
        end

        switch meas(k)
            case 1
                H0=[1,0;0,0];   % only u measured
            case 2
                H0=[0,0;0,1];   % only v measured
            case 3
                H0=eye(2);      % both measured
            otherwise
                H0=zeros(2);    % no measurement if passive, completely failed or
                % not hit
        end
        % the directions of the coordinates u and v are depentent on the radial
        % direction of the intersection point

        % Step 1: derive polar angle Phi
        Phi=atan(paramf(2)/paramf(1));
        if paramf(1)<0
            Phi=pi+Phi;
        elseif paramf(2)<0
            Phi=2*pi+Phi;
        end

        % Step 2: add angles defining the directions
        Psi1=Phi+deltau(fnow);
        Psi2=Phi+deltav(fnow);
        % make sure that both are between 0 and 2*pi
        while Psi1>2*pi
            Psi1=Psi1-2*pi;
        end
        while Psi1<0
            Psi1=Psi1+2*pi;
        end
        while Psi2>2*pi
            Psi2=Psi2-2*pi;
        end
        while Psi2<0
            Psi2=Psi2+2*pi;
        end

        % H matrix of current layer, dependent on intersection point
        Hk=[cos(Psi1) sin(Psi1) 0 0 0
            cos(Psi2) sin(Psi2) 0 0 0];

        paramh=Hk*paramf'; % transformation vector space -> measurement space
        u=paramh(1);
        v=paramh(2);

        if fdistr(fnow)==0    % simulate gaussian distributed errors

            sigsmear(1)=sigu(fnow);
            sigsmear(2)=sigv(fnow);
            % random measurement errors according to computed sigmas
            smear(1)=randn(1)*sigsmear(1);
            smear(2)=randn(1)*sigsmear(2);
            meas_true=[u;v];
            meas_store=H0*(meas_true+smear'); % errors are added
            % H0 eliminates inefficient measurements by neglecting the
            % measurement error. So the 'measurement' of an inefficient
            % layer is the true intersection point. However, in the
            % reconstruction the measurements marked as failed do
            % not contribute!
            coord1(k)=meas_store(1);   % stores u,v
            coord2(k)=meas_store(2);
            
        elseif fdistr(fnow)==1    % simulate uniformly distributed errors

            smear(1)=-du(fnow)/2+rand(1)*du(fnow);  % random measurement errors
            smear(2)=-dv(fnow)/2+rand(1)*dv(fnow);
            meas_true=[u;v];
            meas_store=H0*(meas_true+smear'); % errors are added
            % concerning H0 see comment above
            coord1(k)=meas_store(1);   % stores u,v
            coord2(k)=meas_store(2);
            sigsmear=[du(fnow) dv(fnow)]/sqrt(12);
            
        else
            
            sigsmear=[NaN,NaN];

        end        % if fdistr(k)==1 - end
    end
    
    V(k,:,:)=diag(sigsmear.^2);
    H(k,:,:)=Hk;
    parami=paramf;    

end