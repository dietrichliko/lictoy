function [rms,hist]=mcrms_test(Radius,MC_res,res_true,param_start,param_fit,MCpullhit)

% function rms
% Called by LDT_main
% Main program: LIC_Detector_Toy
%
% Input:    MC_res      Array of residuals at the inner side of the
%                       beamtube, for every track and every coordinate
%           res_true    Residuals between the true simulated parameters and
%                       the fitted parameters at the inner side of the
%                       beamtube, for every track and coordinate
%                       (Phi,z,theta,beta,kappa)
%           param_start Simulated start parameters at inner side of
%                       beamtube, for every track
%           Radius      Radius of beamtube
%           MCpullhit   Logical array, which indicates the tracks for those
%                       Monte-Carlo pulls can be computed
% Output:   hist        arrays for histogramming MC-pulls
%           rms         rms values of hist
%                       rms(1) => RPhi
%                       rms(2) => z
%                       rms(3) => theta
%                       rms(4) => phi
%                       rms(5) => dpt/pt
%                       rms(5) => dpt/pt^2
%
%   RMS calculates the pull quantities at the inner side of the
%   beamtube and and delta p_t/(p_t) and delta p_t/(p_t)^2

global convf unit
global Bz
global Flags

L=MCpullhit;
% logical array to pick out the tracks to display and 
% to fill into the histogram

rms=NaN(1,6);

% delta RPhi
RPhi=res_true(:,1)*Radius;
beta=res_true(:,4);

% to prevent 2*pi errors
Lpi=beta>pi/2;
if sum(Lpi)~=0
    beta(Lpi)=beta(Lpi)-2*pi;
end
Lpi=beta<-pi/2;
if sum(Lpi)~=0
    beta(Lpi)=beta(Lpi)+2*pi;
end
 
% compute residual of phi from beta and Phi
phi=beta+res_true(:,1);

%ptcorr(L)=pt(L).*sqrt(sin(param_start(L,3)));

pt(L)=convf*Bz./param_start(L,5);                   % transverse momentum from curvature
pt_fit(L)=convf*Bz./param_fit(L,5);
dpt(L)=(pt_fit(L)-pt(L))./pt(L);
p(L)=pt(L)./sin(param_start(L,3))';
p_fit(L)=pt_fit(L)./sin(param_fit(L,3))';
dpt2(L)=(p_fit(L)-p(L))./p(L);
%dpt2(L)=dpt(L)./pt(L);
%dpt(L)=MC_res(L,5)./param_start(L,5);               % relative deviation dpt/pt
%dpt2(L)=MC_res(L,5)./(param_start(L,5).*pt(L)');    % relative deviation dpt/pt^2


% computing arrays for MC-pulls and rms of MC-Pulls
if sum(L)~=0
%        [RPhi,z,theta,phi,p_t/(p_t),p_t/(p_t)^2]
    hist=[RPhi(L),res_true(L,2),res_true(L,3),phi(L),dpt(L)',dpt2(L)'];
    rms=sqrt(mean(hist.^2,1));
end

