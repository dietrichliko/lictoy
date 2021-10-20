function [SPR,N,varargout]=paramconversion

% function LDT_ParamConversion
% Called by LDT_ReadParameters
% Main function: LDT_main
%
% Input:    Values (contains read out values from the parameter file)
% Output:   classes SPR,N (contain the Start Parameter Range (SPR),
%               Number of Tracks (N)
%
% LDT_ParamConversion fills the values read from the simulation parameter
% sheet into the classes used by the rest of the program.

warning off
global Flags unit Mass octave

if octave
  namep='simulation_parameters_Octave.txt';
else
  namep='simulation_parameters_Matlab.txt';
end
    
if octave % Octave compatible file readout
  fid=fopen(namep);
  d=0;
  while 1
    d=d+1;
    str=fgetl(fid);
    if strcmp(num2str(str),'-1')
        break;
    end
    A=split(str,':');
    S=size(A);
    if S(1)==2;
        Values{d}=A(2,:);
    else
        Values{d}=[];
    end
  end
  fclose(fid);

  k=1;
  barrelnames=Values{k};
  %if barrelnames(end)~=','
  %  barrelnames(end+1)=',';
  %end
  barrelfile=split(barrelnames,',');
  S=size(barrelfile);
  barrelfile=mat2cell(barrelfile,ones(1,S(1)),S(2));
  for d=1:length(barrelfile)
    str=barrelfile{d};
    str(str==' ')=[];
    barrelfile{d}=str;
  end
  %keyboard
  varargout{1}=barrelfile;
  k=k+1;

  %forwardfile=strread(Values{k},'%s','delimiter',',');
  forwardnames=Values{k};
  %if forwardnames(end)~=','
  %  forwardnames(end+1)=',';
  %end
  forwardfile=split(forwardnames,',');
  S=size(forwardfile);
  forwardfile=mat2cell(forwardfile,ones(1,S(1)),S(2));
  for d=1:length(forwardfile)
    str=forwardfile{d};
    str(str==' ')=[];
    forwardfile{d}=str;
  end
  varargout{2}=forwardfile;
  k=k+1;

else % Matlab compatible
  [Text,Values]=textread(namep,'%s %s','delimiter',':');
  k=1;
end

% Solenid magnetic field
%Bz=str2num(Values{k});
%k=k+1;
% mass of the particles
Mass=str2num(Values{k});
k=k+1;
% Number of events
N.Event=str2num(Values{k});
k=k+1;
% Number of tracks per event
N.Track=str2num(Values{k});
k=k+1;
% Run number
%runnumber=str2num(Values{k});
%k=k+1;
% Minimum measurements
%minmeas=str2num(Values{k});
%k=k+1;
% Vertex filename
%vertexfilename=Values{k};
%k=k+1;
% Momentum filename
%momentumfilename=Values{k};
%k=k+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Converting values of start parameter range (SPR) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=k+1;
% Transverse momentum
SPR.Ptmin=str2num(Values{k});
k=k+1;
SPR.Ptmax=str2num(Values{k});
if length(SPR.Ptmin)~=length(SPR.Ptmax)
    error(['Start parameter range: Momentum (Line ',num2str(k),'): Number of minimum values has to equal number of maximum values!']);
end
% SPR.P_abs if Flags.AbsMomentum is set
SPR.P_absmin=SPR.Ptmin;
SPR.P_absmax=SPR.Ptmax;
k=k+1;

% Angular range in theta
SPR.Thetamin=str2num(Values{k});
k=k+1;
SPR.Thetamax=str2num(Values{k});
if length(SPR.Thetamin)~=length(SPR.Thetamax)
    error(['Start parameter range: Theta (Line ',num2str(k),'): Number of minimum values has to equal number of maximum values!']);
end
k=k+1;

% Angular range in beta
SPR.phimin=str2num(Values{k});
k=k+1;
SPR.phimax=str2num(Values{k});
if length(SPR.phimin)~=length(SPR.phimax)
    error(['Start parameter range: phi (Line ',num2str(k),'): Number of minimum values has to equal number of maximum values!']);
end
k=k+1;

%
% Range in x
%SPR.x=sort(str2num(Values{k})/unit);
%if SPR.x(1)==SPR.x(2)
%    SPR.x=SPR.x+[-eps eps];
%end
%k=k+1;
% Range in y
%SPR.y=sort(str2num(Values{k})/unit);
%if SPR.y(1)==SPR.y(2)
%    SPR.y=SPR.y+[-eps eps];
%end
%k=k+1;
% Range in z
%SPR.z=sort(str2num(Values{k})/unit);
%if SPR.z(1)==SPR.z(2)
%    SPR.z=SPR.z+[-eps eps];
%end
%k=k+1;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Converting values of Flags (Flags) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=k+1;
Flags.Simul=1;
Flags.MulSca=0;
Flags.Smear=0;
Flags.Recon=1;
Flags.DispBadTrack=0;
Flags.Chi2=0;
%Flags.DispPull=0;
%Flags.DispRes=0;
%Flags.Vertex=0;
Flags.AbsMomentum=0;
Flags.ThetaSymmetry=0;
%Flags.Momentum=0;
Flags.ScaleDownTPC=0;
%Flags.DispImp=0;

k=k+1;
% Flags to enable symmetry of the theta range
Flags.ThetaSymmetry=str2num(Values{k});
if ((Flags.ThetaSymmetry ~=0 & Flags.ThetaSymmetry ~=1) | length(Flags.ThetaSymmetry)~=1)
    error('Flags ThetaSymmetry: number of input arguments has to be 1, only 0 or 1 are allowed!');
end
k=k+1;

% Flags to enable use of the absolute value of the momentum rather than the
% traverse momentum
Flags.AbsMomentum=str2num(Values{k});
if ((Flags.AbsMomentum ~=0 & Flags.AbsMomentum ~=1) | length(Flags.AbsMomentum)~=1)
    error('Flags AbsMomentum: number of input arguments has to be 1, only 0 or 1 are allowed!');
end
k=k+1;

% Flags to enable simulation
%Flags.Simul=str2num(Values{k});
%if ((Flags.Simul~=0 & Flags.Simul~=1) | length(Flags.Simul)~=1)
%    error('Flags simulation: number of input arguments has to be 1, only 0 or 1 are allowed!');
%end
%k=k+1;
if Flags.Simul==1
    % Flagss concerning simulation, reconstruction and tests will only be
    % set, if simulation is enabled
    
    % Scale down TPC by factor 5
    Flags.ScaleDownTPC=str2num(Values{k});
    if ((Flags.ScaleDownTPC~=0 & Flags.ScaleDownTPC~=1) | length(Flags.ScaleDownTPC)~=1)
        error('Flags Scale down TPC: number of input arguments has to be 1, only 0 or 1 are allowed!');
    end
    k=k+1;
    % Flags to enable multiple scattering
    Flags.MulSca=str2num(Values{k});
    if ((Flags.MulSca~=0 & Flags.MulSca~=1) | length(Flags.MulSca)~=1)
        error('Flags Multiple scattering: number of input arguments has to be 1, only 0 or 1 are allowed!');
    end
    k=k+1;
    % Flags to enable measurement errors
    Flags.Smear=str2num(Values{k});
    if ((Flags.Smear~=0 & Flags.Smear~=1) | length(Flags.Smear)~=1)
        error('Flags Measurement errors: number of input arguments has to be 1, only 0 or 1 are allowed!');
    end
    k=k+1;
    k=k+1;
    % Flags to enable reconstruction
    %Flags.Recon=str2num(Values{k});
    %if ((Flags.Recon~=0 & Flags.Recon~=1) | length(Flags.Recon)~=1)
    %    error('Flags Reconstruction: number of input arguments has to be 1, only 0 or 1 are allowed!');
    %end
    %k=k+1;
    if Flags.Recon==1        % Flags to enable reconstruction
            
            % Flagss concerning reconstruction will only be set if reconstruction is enabled
            % Flags to enable display of bad tracks
            Flags.DispBadTrack=str2num(Values{k});
            if ((Flags.DispBadTrack~=0 & Flags.DispBadTrack~=1) | ...
                    length(Flags.DispBadTrack)~=1)
                error('Flags Display bad tracks: number of input arguments has to be 1, only 0 or 1 are allowed!');
            end
            k=k+1;
            % Flags to enable chi2
            Flags.Chi2=str2num(Values{k});
            if ((Flags.Chi2~=0 & Flags.Chi2~=1) | length(Flags.Chi2)~=1)
                error('Flags Chi2: number of input arguments has to be 1, only 0 or 1 are allowed!');
            end
            k=k+1;
            k=k+1;
            %
            % Flags to enable histograms of pulls
            %Flags.DispPull=str2num(Values{k});
            %if ((Flags.DispPull~=0 & Flags.DispPull~=1) | ...
            %        length(Flags.DispPull)~=1)
            %    error('Flags Pulls histograms: number of input arguments has to be 1, only 0 or 1 are allowed!');
            %end
            %k=k+1;
            % Flags to enable histograms of residuals
            %Flags.DispRes=str2num(Values{k});
            %if ((Flags.DispRes~=0 & Flags.DispRes~=1) | ...
            %        length(Flags.DispRes)~=1)
            %    error('Flags Residuals histograms: number of input arguments has to be 1, only 0 or 1 are allowed!');
            %end
            %k=k+1;
            
            % Flags to enable histograms of impact parameters
            %Flags.DispImp=str2num(Values{k});
            %if ((Flags.DispImp~=0 & Flags.DispImp~=1) | ...
            %        length(Flags.DispImp)~=1)
            %    error('Flags Impact histograms: number of input arguments has to be 1, only 0 or 1 are allowed!');
            %end
            %k=k+1;

            %Flags.Vertex=str2num(Values{k});
            %if ((Flags.Vertex~=0 & Flags.Vertex~=1) | length(Flags.Vertex)~=1)
            %error('Flags Vertex file: number of input arguments has to be 1, only 0 or 1 are allowed!');
            %end
            %k=k+1;
            
            %Flags.Momentum=str2num(Values{k});
            %if ((Flags.Momentum~=0 & Flags.Momentum~=1) | length(Flags.Momentum)~=1)
            %error('Flags Momentum file: number of input arguments has to be 1, only 0 or 1 are allowed!');
            %end
            %k=k+1;
            %}
        
    else
        k=k+2;
    end     % if Flags.Recon==1
end         % if Flags.Simul==1

