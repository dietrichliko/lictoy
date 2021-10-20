function [param_fit,Cov,chi2]=testreconstruction(paramr,Hk,Vk,Ak,coord1,...
    coord2,meas,reftype,refindex,varMS)

% function TESTRECONSTRUCTION
% Called by LTD_main
% Main program: LTD_main
%
% Input:    paramt      Start parameters at inner side of beamtube
%           Hk_store    Array of H matrices for every layer
%           Vk_store    Array of measurement error matrices
%           coord1      RPhi or u measurements of each layer, respectively
%           coord2      z or v measurements of each layer, respectively
%           eff         Indicator for inefficiencies (0: no measurement, 1:
%                       only RPhi, 2: only z, 3: both) of each layer
%           hit         Indicator which layer was hit
%           cyl         Indicates which of the layers are cylinder layers
%           mass        Mass of the particle
%
% Output:   param_fit   Fitted state vectors of each layer
%           Cf_store    Corresponding covariance matrices of the fitted
%                       state vectors
%           chi2        Chi2 of the fit
%
% TESTRECONSTRUCTION performs the following steps for a single track: 
% First the reference track (the expansion points for the Kalman filter) 
% is computed.
% The reference track is the completely undisturbed extrapolation of the
% start parameters through the whole detector, without multiple scattering.
% Having derived the reference track and the according derivative matrices
% the Kalman filter is performed with parameters, which are the local 
% difference between the reference track and the real measurements (linear 
% model between two detectors). The filter gives back the fitted parameters
% and corresponding covariance matrices for every layer, which are the 
% output parameters of the function IRECONSTRUCTION.
% The chi2 of the filter is only computed if desired, because it is the
% most time consuming part of the filter.
%

global Flags unit ITrack IEvent fidlog convf Bz radius hitpattern
warning off

%------------------------ Reconstruction ---------------------------------

chi2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------- Kalman filter -----------------------------

Ce=eye(5);              % initial covariance matrix and initial
parame=zeros(1,5);      % parameters to get the filter going

measNr1=size(paramr);
measNr1=measNr1(1);
%
if reftype(measNr1)
    % barrel initialization
    Ce(1,1)=(20/unit/radius.b(refindex(measNr1)))^2;
    Ce(2,2)=(20/unit)^2;
    Ce(3,3)=(2*1e-2)^2;
    Ce(4,4)=(20*1e-2)^2/4;
    Ce(5,5)=(convf/2)^2;
    %Ce(5,5)=(1e-4*33)^2;
    Ce=Ce*100;
else
    %}
    pt=convf*Bz/paramr(end,5);
    % forward initialization
    Ce(1,1)=(5/unit)^2;
    Ce(2,2)=(5/unit)^2;
    Ce(3,3)=(20*1e-2)^2;
    Ce(4,4)=(20*1e-2)^2;
    %Ce(5,5)=((paramr(end,5)*pt/0.1)/4)^2;
    Ce(5,5)=(convf*10)^2;
    %Ce(5,5)=(1e-4*33)^2;
    %Ce=Ce*10;
end

% All filtering is done with parameters which refer to the
% the reference track and the real parameters (local linear model).
% Starting outside the outermost considered layer, first multiple
% scattering is added, giving the parameters at the inner side of the
% layer. The measurements are assumed to be located at the inner side.
% Here the famous 'weighted mean' of the Kalman-filter is calculated.
% The last step is the propagation to the next layer.
% If desired, the total chi2 of the fit is computed, too.

%pr=abs((convf*Bz)/(paramr(1,5)*sin(paramr(1,3))));
pr=abs((convf*Bz)/paramr(1,5)/sin(paramr(1,3)));
% Momentum of reference track

for k=measNr1:-1:1    % filtering backwards until first layer hit
    
    % Add multiple scattering, blowing up the covariance matrix
    if Flags.MulSca
        fac=-convf*Bz/pr*cos(paramr(k,3))/sin(paramr(k,3))^2;
        Ce(3,3)=Ce(3,3)+varMS(k);
        Ce(4,4)=Ce(4,4)+varMS(k)/sin(paramr(k,3))^2;
        Ce(3,5)=Ce(3,5)+fac*varMS(k);
        Ce(5,3)=Ce(3,5);
        Ce(5,5)=Ce(5,5)+fac^2*varMS(k);
    end
    
    if meas(k)   % whole filtering
        % procedure only for existing
        % measurements and traversed layers

        V=squeeze(Vk(k,:,:));% meas. error matrix of current layer
        H=squeeze(Hk(k,:,:));% H matrix of current layer
        
        switch meas(k)
            case 1
                H(2,:)=0; H0=[1,0;0,0];  % only coord1 measured
            case 2
                H(1,:)=0; H0=[0,0;0,1];  % only coord2 measured
            otherwise
                H0=eye(2); % no changes, both coordinates measured
        end
		
        if reftype(k)
            % to prevent 2*pi-errors
            if (coord1(k)/radius.b(refindex(k))-paramr(k,1))>pi
                paramr(k,1)=paramr(k,1)+2*pi;
            elseif (coord1(k)/radius.b(refindex(k))-paramr(k,1))<-pi
                paramr(k,1)=paramr(k,1)-2*pi;
            end
        end

        paramh=H*paramr(k,:)';
        % transforms reference parameters to measurement space
        mk(1)=coord1(k)-paramh(1);   % difference between reference track and
        mk(2)=coord2(k)-paramh(2);   % measurements is used in the filtering 
                                % procedure
       
        %mk=H0*mk';
        %mk=mk';
                                
        % Kalman filter equations
        Ek=V+H*Ce*H';
        Kk=Ce*H'*inv(Ek);                      % gain matrix
        paramf=(parame'+Kk*(mk'-H*parame'))';  
        % fitted parameters (weighted mean)
        Cf=(eye(5)-Kk*H)*Ce;                   
        % cov. matrix of fitted parameters
        
        if Flags.Chi2==1 

            mk(1)=coord1(k);             % u measurement of current layer
            mk(2)=coord2(k);             % v measurement of current layer
            paramh=paramf+paramr(k,:);  % reference parameters added -> 
                                        % real parameters
                                                                             
            switch meas(k)
                case 1
                    H0=[1,0];   % only RPhi measured
                case 2
                    H0=[0,1];   % only z measured
                otherwise
                    H0=eye(2);      % both measured
            end
                                        
            % With help of H0 inefficient measurements are neglected, so
            % that they don't contribute to the total chi2

            chi2res=H0*mk'-H0*(H*paramh'); 
            % residual, inefficient measurements terminated by H0
            Ek=H*Cf*H';           % vector space -> measurement space
            Rk=H0*V*H0'-H0*Ek*H0'; % matrix for normalization
            chi2plus=chi2res'*(Rk^-1)*chi2res;  % increase of chi2
            chi2=chi2+chi2plus;   % increasing of total chi2
        end     % if Flags.Chi2==1 - end

    else        % for totally inefficient and passive layers
        paramf=parame;        % no weighted mean in a passive scatterer
        Cf=Ce;                % no update of the covariance matrix
    end         % if meas(k) - end

    Cov(k,:,:)=Cf;               % stores the covariance matrix
    param_fit(k,:)=paramf+paramr(k,:);% stores the real fitted parameters

    if k==1
        break               % last propagation skipped when reached
    end                     % the inner side of the beam tube

    warning off
	Fk=inv(squeeze(Ak(k,:,:)));   
    warning on
	parame=paramf*Fk';      
	Ce=Fk*Cf*Fk';           
         
end    % loop over the layers - end
%keyboard