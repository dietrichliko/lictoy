function [RPhi_store,z_store,Vk_store,eff_ind]=tpc_ineff(hit,paramh,paramf,...
    tpc_zmax,tpc_zmin,eff1,eff2,sigRPhia,sigRPhib,sigza,sigzb)

% Function tcp_ineff
% called by isimualtion 
%
%   input:  paramf			parameters of the reference track
%			paramh			parameters of the reference track in the
%							measurment space
%			tpc_zmin/zmax	length coordinates of the TPC
%           eff1			efficiency in RPHI
%           eff2			efficiency in z 
%			sigRPhia/b		 
%			sigza/b		simgas in RPhi and z
%
%	output:	eff_ind			indicates which coordinate was measured
%			RPhi_store		measurements in RPhi
%			z_store			measurements in z
%			Vk_stroe		error matrix
%			
%----------------------------------------------------------------------

	global unit

    meas1=rand(1)<=eff1;
    
    if eff2==-1
        meas2=meas1;
    else
        meas2=rand(1)<eff2;
    end

    if meas1 & meas2
            eff_ind=3;
    elseif meas2
            eff_ind=2;
    elseif meas1
            eff_ind=1;
    else
            eff_ind=0;
    end

    % only measurements when layer was hit 
    switch eff_ind
            case 1
                H0=[1,0;0,0];   % only RPhi measured
            case 2
                H0=[0,0;0,1];   % only z measured
            case 3
                H0=eye(2);      % both measured
			otherwise
				H0=zeros(2);	% no measurement
	end	
	if hit==0
        H0=zeros(2);
    end	
	%----------------------------------------------------------------------
        
	% simulate gaussian distributed errors
    % formula for errors in TPC is also usable for other detector layers
    % if sigRPhib=sigzb=0
    % sigma^2 = [ sigma1^2 + sigma2^2 * (abs( z - zmax )) ]
    if paramf(2)>0
		zdiff=abs(paramf(2)-tpc_zmax);
	else
		zdiff=abs(paramf(2)-tpc_zmin);
	end
    
	% zdiff in [cm] for the diffusion formula
	zdiff=zdiff*1e-1*unit; 
    sigsmear(1)=sqrt((sigRPhia*1e-3*unit)^2+((sigRPhib*1e-3*unit)^2)*zdiff);
    sigsmear(2)=sqrt((sigza*1e-3*unit)^2+((sigzb*1e-3*unit)^2)*zdiff);
    
	% sigsmear in [m]
    sigsmear=sigsmear*1e3/unit;  % in [mm]
   
	
	% random measurement errors according to computed sigmas
    smear=randn(1,2).*sigsmear;
    meas_true=[paramh(1);paramh(2)];
    meas_store=meas_true+H0*smear'; % errors are added
	
    RPhi_store = meas_store(1);         % stores RPhi and z
    z_store    = meas_store(2);
    Vk_store   = diag(sigsmear.^2);     % stores the error matrix