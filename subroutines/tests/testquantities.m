function [MC_pull,MC_res,pull,res,res_true,MC_chi2]=testquantities(paramt,...
    paramf,coord1,coord2,Vk,Hk,Cf_store,meas,Radius)

% function BTESTQUANTITIES
% Called by LDT_main
% Main program: LDT_main
%
% Input:    paramt      True parameters at inner side of beamtube and at
%                       first detector layer
%           paramf      Fitted parameters at inner side of beamtube and at
%                       first detector layer
%           coord1        coord1 measurement at first detector layer
%           coord2           coord2 measurement at first detector layer
%           Vk          Measurement error matrices at inner side of
%                       beamtube and at first detector layer
%           Hk          H matrices at inner side of beamtube and at
%                       first detector layer
%           Cf_store    Covariance matrices of the fitted parameters at
%                       inner side of beamtube and at first detector layer
%           pullhit     Indicator if pulls at first detector layer can be
%                       computed, or if the layer was inefficient or not
%                       hit for the present track
%           Radius      Radii of beamtube and first detector layer
%
% Output:   MC_pull     Monte Carlo pull quantities, computed from start
%                       parameters and the fitted parameters at the inner
%                       side of the beamtube, for every coordinate
%                       (Phi,coord2,theta,beta,kappa)
%           MC_res      Monte Carlo residuals at the inner side of the
%                       beamtube, for every coordinate
%           pull        Real pull quantities at innermost detector layer,
%                       computed from the fitted parameters and the real
%                       measurements, for every coord. (coord1,coord2')
%           res         Real residuals at the innermost detector layer, for
%                       every coordinate
%           res_true    Residuals between the parameters of the true
%                       intersection points and the fitted parameters, for
%                       every coordinate (Phi,coord2,theta,beta,kappa)
%
% BTESTQUANTITIES computes the test quantities of the barrel fit. At the
% inner side of the beamtube the Monte-Carlo (MC) residuals and pulls are
% computed using the simulated start vectors and the fitted parameters. At
% the first detector layer the real residuals and pulls are computed using
% the real (simulated) measurements and the fitted measurements.
%

%warning off
global limpull Flags ITrack radius
warning off;

pull=zeros(1,2);    % real pull quantities at innermost detector
res=zeros(1,2);     % real residuals at innermost detector
res_true=zeros(1,5);% residuals (param_true-param_fit) for each coordinate


% Monte Carlo pulls
Cf=squeeze(Cf_store(1,:,:)); % covariance matrix of current track at inner
% side of beamtube
% to prevent 2*pi-errors
% for Phi
if (paramf(1,1)-paramt(1,1))>pi
    paramf(1,1)=paramf(1,1)-2*pi;
elseif (paramf(1,1)-paramt(1,1))<-pi
    paramf(1,1)=paramf(1,1)+2*pi;
end
% for beta (bzw. phi)
if (paramf(1,4)-paramt(1,4))>pi
    paramf(1,4)=paramf(1,4)-2*pi;
elseif (paramf(1,4)-paramt(1,4))<-pi
    paramf(1,4)=paramf(1,4)+2*pi;
end
% for theta
if (paramf(1,3)-paramt(1,3))>pi
    paramf(1,3)=paramf(1,3)-2*pi;
elseif (paramf(1,3)-paramt(1,3))<-pi
    paramf(1,3)=paramf(1,3)+2*pi;
end
%}

% Monte-Carlo-pulls compare the fitted parameters
% inside the beam tube with the real simulated parameters at BT
MC_res(1:5)=paramf(1,:)-paramt(1,:);    % computing residuals
MC_pull(1:5)=MC_res./sqrt(diag(Cf))';   % computing pulls
warning off;
MC_chi2=MC_res*inv(Cf)*MC_res';         % Monte Carlo chi^2
%warning on

MC_res(1)=Radius(1)*MC_res(1);          % residuals dim=[length]


% Real pulls
% Real residuals and pulls at first detector, using the fitted
% "measurements" and the simulated measurements
if meas
    Cf=squeeze(Cf_store(2,:,:)); % covariance matrix of current track at
    % innermost detector layer

    Rk=Vk-Hk*Cf*Hk';          % matrix for normalization of pulls
    paramh=Hk*paramf(2,:)';   % transf. state vector space -> meas. space

    if length(Radius)>1
        % to prevent 2*pi-errors for RPhi
        if (paramf(2,1)-coord1/Radius(2))>pi
            paramf(2,1)=paramf(2,1)-2*pi;
        elseif (paramf(2,1)-coord1/Radius(2))<-pi
            paramf(2,1)=paramf(2,1)+2*pi;
        end
    end

    res(1)=paramh(1)-coord1;			% residual of coord1
    res(2)=paramh(2)-coord2;				% residual of coord2

    vars=diag(Rk);            % variances for normalization
    pull=res./sqrt(vars)';    % pulls
    
    if meas==1
        pull(2)=nan;
    elseif meas==2
        pull(1)=nan;
    end
else %
    pull(:)=NaN;
end

% Residuals between true tracks (incl. MS) and fitted tracks
res_true=paramt(1,:)-paramf(1,:);