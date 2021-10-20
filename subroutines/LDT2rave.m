function LDT2rave(filename,varargin)
%(SPR,param_start,param_fit,R,Cf_store)
% function rave
% Called by output
% Main program: LDT_main
%
% Input:    SPR         Structure containing the start parameter range
%           param_start Generated start parameters at inner side of
%                       beamtube, for every fitted track
%           param_fit   Fitted parameters at inner side of beamtube, for
%                       every fitted track
%           R           Radius of beamtube
%           Cf_store    Covariance matrices at inner side of beamtube, for
%                       every fitted track
%
% Output:   none
%
% RAVE delivers an output of the fitted parameters at the inner side
% of the beamtube in a text file. The parameters are converted to cartesian
% using the function VERTCONV and are then written to a text file in the
% DATA HARVESTER's comma seperated variables format. This textfile can be
% read out with said data harvester and fed into the RAVE/VERTIGO vertex
% fit toolkit.
% To monitor the correctness of the conversion, the pulls and the chi2
% after the conversion are computed, if desired.

%global Bz unit Flags
%global fidlog hwait

warning off;
if isempty(filename)
    error('Please enter a filename!')
end

if isempty(varargin)
    fnameresults='results.mat';
else
    fnameresults=varargin{1};
end

disp(['Using results file: ',fnameresults]);
load(fnameresults,'vertex','param_start','param_fit','R','Cf_store','Bz','unit','Flags','convf');

[N,d]=size(param_start);
param6=zeros(N,6);		% cartesian fitted parameters for every track
MCparam6=zeros(N,6);	% cartesian start parameters for every track
C6=zeros(N,6,6);		% 6x6 covariance matrix, for every track
pull6=zeros(N,6);		% pulls after conversion

%fidlog=fopen('LIC_Detector_Toy.log','a');
disp('Conversion to cartesian coordinates (x,y,z,p_x,p_y,p_z)');

for i=1:N
    % conversion of fitted parameters and covariance matrix (at BT)
    parami=squeeze(param_fit(i,:));
    Cp=squeeze(Cf_store(i,:,:));
    % vertconv does conversion of parameters and covariance matrix
    [qaramfit,Cq]=vertconv(Bz,R,parami,Cp,abs(convf),unit); 
    param6(i,:)=qaramfit;
    C6(i,:,:)=Cq;
    % conversion of simulated Monte-Carlo parameters (at BT)
    parami=squeeze(param_start(i,:));
    [qaramMC,Dummy]=vertconv(Bz,R,parami,eye(5),abs(convf),unit);
    MCparam6(i,:)=qaramMC;
	
	
    if Flags.Recon
        if Flags.Chi2
			% chi2 only for the momenta
            V=squeeze(Cq(4:6,4:6));
			chi2res=qaramfit(4:6)-qaramMC(4:6);		% chi2 residual
            chi2(i)=chi2res*inv(V)*chi2res';		% chi2 test quantity
		end

        res=qaramfit-qaramMC;						% residuals
        pull6(i,:)=res./sqrt(diag(Cq))';			% pull quantities
	end
end

% output of test quantities of conversion
if Flags.Recon
    meanpull6=mean(pull6);			% mean of pulls
    stdpull6=std(pull6);			% standard deviation of pulls	
	
    disp('Test quantities after conversion 5dim -> 6dim');
    if Flags.Chi2
        meanChi2_6=mean(chi2);
        disp('Mean of chi2 (only momentum coordinates considered, should be 3)');
        printf('%8f\n',meanChi2_6);
    end
    disp('Pulls of  x           y           z          px          py          pz');
    disp(['mean: ',num2str(meanpull6)]);
    disp(['std:  ',num2str(stdpull6)]);
    
	
end

% Output values for vertex fit
VTX=[vertex(1)/10*unit,vertex(2)/10*unit,vertex(3)/10*unit];					% vertex position in [cm]
Q=-sign(param_fit(:,5)*Bz);						% charge of particles
vertprint(N,Q,VTX,MCparam6,param6,C6,filename);			% creates output
%fclose(fidlog);