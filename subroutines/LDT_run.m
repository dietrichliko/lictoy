function LDT_run(varargin)
%
% Close histogram figures if they exist
%save _temp.mat action nameb namef parameterfile messagehandle warninghandle;
%clear;
%load _temp.mat action nameb namef parameterfile messagehandle warninghandle;
%delete _temp.mat

%path(path,'subroutines');
%path(path,'subroutines/gui');
%path(path,'subroutines/outputs');
%path(path,'subroutines/simulation');
%path(path,'subroutines/tests');

%switch action
%    case 'display2D'
%        display=1;
%    case 'display3D'
%        display=2;
%    otherwise
%        display=0;
%end

tic;
warning off;
%warning on
format short;
format loose;

%offset=10e8;
offset=1e5;
% to reproduce the same sequence of random numbers.
rand('state',offset);
randn('state',offset);

global Flags GeomVersion  % GeomVersion for version control of *.geom files
global Bz Mass convf sinbmx SPR limpull unit RadiusIR ITrack IEvent
global hw_event hw_track numhitnr1 hwait TPC
global Flags fidlog N radius zpos delta Xlen detnr1
global eff distr sig d delta name hitpattern octave 
if ~octave
    global whandle mhandle interrupt
    % setting up default properties for graphical output
    set(0,'defaultfigureunits','pixels')
    set(0,'defaulttextfontsize',10)
    set(0,'defaultaxesfontsize',9)
    %set(0,'defaulttextfontname','times')
    %set(0,'defaultaxesfontname','times')
end
interrupt=0;

%------------ Close histogram windows ---------------------------------
%disfig=1;
%hisfig=2:4;
%if ~display close(hisfig(ishandle(hisfig))); end

%------------  Clear message and warning handles ----------------------
%whandle=warninghandle;
%mhandle=messagehandle;

%------------ Open log file -------------------------------------------
fidlog=fopen('LDT.log','w');
%str=['Action: ',action];
%fprintf(fidlog,'%s\n',str);

%unit=1;                        % unit=1   : computation in [mm]
% unit=10  : computation in [cm]
% unit=1000: computation in [m]
RadiusIR=1e-3/unit;
convf=0.299792458*1e-3*unit;    % conversion factor [Gev/c T^(-1) m^(-1)]
sinbmx=0.9;                     % maximal value of sin(beta) allowed in propagation
pt_limit=0.5;                   % lower limit of the transverse momentum
limpull=3;                      % pulls > limpull are considered bad
numhitnr1=0;

%--------------------------------------------------------------------------
%
%[SPR,N]=LDT_ReadParameters(parameterfile);
%[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=LDT_ReadGeometry(barrelfile,forwardfile);
%merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2);
%}
minmeas=3;

% program version
PV=1;
% run number
runnumber=2020;
% vector which contains numbers to set state of the random generators
%rand_init=[PV,GeomVersion.BNo*100,GeomVersion.FNo*100,runnumber];
if octave
    [SPR,N,nameb,namef]=paramconversion;
else
    [SPR,N]=paramconversion;
    nameb=varargin{1};
    namef=varargin{2};
    mhandle=varargin{3};
    whandle=varargin{4};
    %------------  Clear message and warning handles ----------------------
    str=get(mhandle,'String');
    str{1}='Messages:';
    str{2}='No messages';
    str{3}=' ';
    str{4}=' ';
    set(mhandle,'String',str);
    str=get(whandle,'String');
    str{2}='No warnings';
    str{3}=' ';
    str{4}=' ';
    set(whandle,'String',str);
    set(whandle,'foregroundcolor',[0 0.5 0])
    drawnow;
end

%
%if display
%    %[SPR,N]=LDT_ReadParameters(parameterfile);
%    [VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=LDT_ReadGeometry(nameb{1},namef{1});
%    merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2);
%    [theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
% angles to seperate forward,
% intermediate, barrel and rear
%    varargout{1}=geometry3D(nameb{1},namef{1},VTX,SIT,TPC,SET,FM1,...
%        FM2,RM1,RM2,SPR,theta1,theta2,theta3,theta4,display);
% calls function 'geometry', if desired,
% which draws a sketch of the arrangement
%    fclose(fidlog);
%    return
%end
%}


if Flags.Simul==1
    %tic
    N.Pt=length(SPR.Ptmin);
    N.Theta=length(SPR.Thetamin);
    N.barrel=length(nameb);
    N.forward=length(namef);
    if N.barrel~=N.forward
        error('Please select the same number of barrel geometries and forward geometries');
    end

    if N.barrel>1
        N.Curve=N.barrel;
        if N.Pt>=N.Theta
            % compute in dependence on pt
            Flags.Curves=5;
            N.Point=N.Pt;
        else
            % compute in dependence on theta
            Flags.Curves=6;
            N.Point=N.Theta;
        end
    else
        if ( N.Pt>1 & N.Theta>1 )
            N.Point=max(N.Pt,N.Theta);
            N.Curve=min(N.Pt,N.Theta);
            if N.Pt<N.Theta
                Flags.Curves=4; % calculate curves in dependence on theta
            else
                Flags.Curves=3; % calculate curves in dependence on pt
            end
        elseif ( N.Pt>1 | N.Theta>1 )
            N.Point=max(N.Pt,N.Theta);
            N.Curve=1;
            Flags.Curves=1;
            if max(N.Pt,N.Theta)==N.Pt
                Flags.Curves=2;
            else
                Flags.Curves=1;
            end
        else
            N.Point=1;
            N.Curve=1;
            Flags.Curves=0;
        end
    end

    barrelfile=nameb{1};
    forwardfile=namef{1};
    [VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,nameb{1},namef{1}]=geomconversion(barrelfile,forwardfile);
    mcount=merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,barrelfile,forwardfile);
    [theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
    % angles to seperate forward,
    % intermediate, barrel and rear
    
    if octave
        count=1;
        printf('\n%s\n','Progress...'); fflush(stdout);
    end
    for ICurve=1:N.Curve

        rand('state',offset);    
        randn('state',offset);
        
        if (Flags.Curves==5 | Flags.Curves==6) & ICurve>1
            barrelfile=nameb{ICurve};
            forwardfile=namef{ICurve};
            [VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,nameb{ICurve},namef{ICurve}]=geomconversion(barrelfile,forwardfile);
            mcount=merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,barrelfile,forwardfile);
            [theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
            % angles to seperate forward,
            % intermediate, barrel and rear
        else
            barrelfile=nameb{1};
            forwardfile=namef{1};
        end
        
        %[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=geomconversion(barrelfile,forwardfile);
        %mcount=merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,barrelfile,forwardfile);
        %[theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
        % angles to seperate forward,
        % intermediate, barrel and rear
        
        start.x=sort(SPR.x);
        start.y=sort(SPR.y);
        start.z=sort(SPR.z);
        start.phi=sort([SPR.phimin,SPR.phimax]);
        
        for IPoint=1:N.Point
            
            if ~octave
            rand('state',offset);    
            randn('state',offset);
            
            str=get(mhandle,'String');
            if Flags.Curves
                str{2}=['Curve ',num2str(ICurve),' of ',num2str(N.Curve),', point ',num2str(IPoint),' of ',num2str(N.Point)];
            else
                str(2)=[];
            end
            set(mhandle,'String',str);drawnow;
        end
        
            
            switch Flags.Curves
                case 6
                    start.Theta(1)=SPR.Thetamin(IPoint);
                    start.Theta(2)=SPR.Thetamax(IPoint);
                    start.Pt(1)=SPR.Ptmin(1);
                    start.Pt(2)=SPR.Ptmax(1);
                case 5
                    start.Theta(1)=SPR.Thetamin(1);
                    start.Theta(2)=SPR.Thetamax(1);
                    start.Pt(1)=SPR.Ptmin(IPoint);
                    start.Pt(2)=SPR.Ptmax(IPoint);
                case 4
                    start.Theta(1)=SPR.Thetamin(IPoint);
                    start.Theta(2)=SPR.Thetamax(IPoint);
                    start.Pt(1)=SPR.Ptmin(ICurve);
                    start.Pt(2)=SPR.Ptmax(ICurve);
                case 3
                    start.Theta(1)=SPR.Thetamin(ICurve);
                    start.Theta(2)=SPR.Thetamax(ICurve);
                    start.Pt(1)=SPR.Ptmin(IPoint);
                    start.Pt(2)=SPR.Ptmax(IPoint);
                case 2
                    start.Theta(1)=SPR.Thetamin(ICurve);
                    start.Theta(2)=SPR.Thetamax(ICurve);
                    start.Pt(1)=SPR.Ptmin(IPoint);
                    start.Pt(2)=SPR.Ptmax(IPoint);
                case 1
                    start.Theta(1)=SPR.Thetamin(IPoint);
                    start.Theta(2)=SPR.Thetamax(IPoint);
                    start.Pt(1)=SPR.Ptmin(ICurve);
                    start.Pt(2)=SPR.Ptmax(ICurve);
                case 0
                    start.Theta(1)=SPR.Thetamin(ICurve);
                    start.Theta(2)=SPR.Thetamax(ICurve);
                    start.Pt(1)=SPR.Ptmin(IPoint);
                    start.Pt(2)=SPR.Ptmax(IPoint);
            end

            %fidlog=fopen('LIC_Detector_Toy.log','w');
            %str=['Action: ',action];
            %fprintf(fidlog,'%s\n',str);

            % Variables to collect data from different events
            AMCpullhit=[];Apullhit=[];Aparam_start=[];ALbarrel=[];Ameas=[];
            ALforward=[];ALinterm=[];Avertex=[];AMC_tracknr=[];Aoffset=[];
            Aparam_true=[];Acoord1=[];Acoord2=[];Ahitpattern=[];Aregion=[];
            fnum=0;                             % number of forward tracks
            rnum=0;                             % number of rear tracks
            bnum=0;                             % number of barrel tracks
            inum=0;                             % number of intermediate tracks

            if Flags.Recon
                Apulldet=[];Apulltype=[];AMC_pull=[];AMC_res=[];Apull=[];Ares=[];
                Ares_true=[];ACf_store=[];Aparam_fit=[];Achi2=[];Andf=[];
                AMC_chi2=[];Aresout_coordinate=[];Achi_hist=[];Achi22=[];
                %Abad=zeros(1,5); Abadcount=0;
                %Aparam_fit_alldets=[];
                %ACf_store_alldets=[];
                badcount=0;                         % counts total number of bad tracks
                bad=zeros(1,5);                     % counts bad pulls seperately
            end

            %Generate waitbars
            Flags.wait=0;
            if Flags.wait waitbars; end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%% Eventloop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for IEvent=1:N.Event

                % Initialisation of variables
                chi2_hist=zeros(N.Track,2);
                resout_coordinate=zeros(N.Track,1); % to store the fixed coord. of the detector to compute the pulls
                vertex=zeros(3,1);                  % Vertex coordinates
                pullhit=zeros(1,N.Track);           % Array for real pulls
                MC_pull=zeros(N.Track,5);           % Array for MC pulls
                MC_res=zeros(N.Track,5);            % Array for MC residuals
                pull=zeros(N.Track,2);              % Array for real pulls
                res=zeros(N.Track,2);               % Array for real residuals
                res_true=zeros(N.Track,5);          % Array for
                param_start_all=zeros(N.Track,5);   % Array for start parameters
                Cf_store_all=zeros(N.Track,5,5);    % Array for covariance matrices
                chi2=zeros(N.Track,1);              % Array for chi2s of the fit
                pulltype=ones(1,N.Track);           % Array for pull types (1=barrel,0=forward)
                MCpullhit=logical(ones(1,N.Track)); % Logical array for MC pulls
                Lbarrel=logical(zeros(1,N.Track));  % Logical array for barrel tracks
                Lforward=logical(zeros(1,N.Track)); % Logical array for forward tracks
                Linterm=logical(zeros(1,N.Track));  % Logical array for interm. tracks
                param_fit_all=zeros(N.Track,5);     % Array for fitted parameters
                param_true_all=zeros(N.Track,5);    % Array for true (simulated) parameters
                coord1_all=zeros(N.Track,2);        % Array for first simulated coordinate
                coord2_all=zeros(N.Track,2);        % Array for second simulated coordinate
                hitpattern_all=cell(N.Track,2); % Array for hit pattern
                region=zeros(1,N.Track);
                meas_all=zeros(N.Track,N.BLayer+max(N.FLayer,N.RLayer));
                pulldet=zeros(1,N.Track);
                ndf=zeros(1,N.Track);
                MC_chi2=zeros(1,N.Track);
                % counts the total number of bad pulls for every coordinate
                % (Phi,z,theta,beta,1/R).
                % Pulls greater than limpull (usually 3) are considered bad
                %barrel=ones(1,N.Track);             % Indicator for barrel tracks

                % hook for external access event by event
                % enter Vertex (x,y,z) in dependence on IEvent

                param_IR(1)=start.x(1)+rand(1)*(start.x(2)-start.x(1));
                param_IR(2)=start.y(1)+rand(1)*(start.y(2)-start.y(1));
                param_IR(3)=start.z(1)+rand(1)*(start.z(2)-start.z(1));

                vertex=[param_IR(1:3)];

                %--------------------------------------------------------------------------
                %--------------------- Trackloop ------------------------------------------
                %--------------------------------------------------------------------------

                for ITrack=1:N.Track
                    % hook for external access track by track
                    % enter track parameters in dependence on ITrack
                    clear hitpattern;
                    warning off;
                    
                    if octave
                        str=[];
                        if ~mcount
                            for k=1:count
                                str=[str,'\b'];
                            end
                        else
                            str='';
                            mcount=0;
                        end
                        printf(str);
                        fflush(stdout);
                        %fprintf('\r');
                        %clc;
                        if Flags.Curves
                            if N.Event>1
                                count=printf('Curve %d of %d, point %d of %d, event %d of %d, track %d of %d     ',...
                                    ICurve,N.Curve,IPoint,N.Point,IEvent,N.Event,ITrack,N.Track);
                            else
                                count=printf('Curve %d of %d, point %d of %d, track %d of %d     ',...
                                    ICurve,N.Curve,IPoint,N.Point,ITrack,N.Track);
                            end
                        else
                            if N.Event>1
                                count=printf('Event %d of %d, track %d of %d      ',...
                                    IEvent,N.Event,ITrack,N.Track);
                            else
                                count=printf('Track %d of %d     ',...
                                    ITrack,N.Track);
                                %message=['Track ',num2str(ITrack),' of ',num2str(N.Track)];
                            end
                        end
                        %waitbar(ITrack/N.Track)
                        %printf(message);
                        fflush(stdout);
                    else
                        str=get(mhandle,'String');
                        if Flags.Curves k=3;
                        else k=2; end
                        if N.Event>1
                            str{k}=['Event ',num2str(IEvent),' of ',num2str(N.Event),', track ',num2str(ITrack),' of ',num2str(N.Track)];
                            set(mhandle,'String',str);drawnow;
                            if Flags.wait waitbar(IEvent/N.Event,hw_event,str{2}); end
                        else
                            str{k}=['Track ',num2str(ITrack),' of ',num2str(N.Track)];
                            set(mhandle,'String',str);drawnow
                        end
                        if Flags.wait waitbar(ITrack/N.Track,hw_track,str{tpos}); end
                    end
                    %disp(num2str(length(message)));
                    %drawnow;

                    %if Flags.Curves k=3;
                    %else k=2; end
                    %if Flags.Curves
                    %    if N.Event>1
                    %        fprintf(['Curve ',num2str(ICurve),' of ',num2str(N.Curve),', point ',num2str(IPoint),' of ',num2str(N.Point),', ','Event ',num2str(IEvent),' of ',num2str(N.Event),', track ',num2str(ITrack),' of ',num2str(N.Track),'\r'])
                    %    else
                    %        disp(['Curve ',num2str(ICurve),' of ',num2str(N.Curve),', point ',num2str(IPoint),' of ',num2str(N.Point),', ','Track ',num2str(ITrack),' of ',num2str(N.Track)])
                    %    end
                    %else
                    %    if N.Event>1
                    %        disp(['Event ',num2str(IEvent),' of ',num2str(N.Event),', track ',num2str(ITrack),' of ',num2str(N.Track)]);
                    %        if Flags.wait waitbar(IEvent/N.Event,hw_event,str{2}); end
                    %    else
                    %        fprintf(['\r','Track ',num2str(ITrack),' of ',num2str(N.Track)]);
                    %    end
                    %end
                    if Flags.wait waitbar(ITrack/N.Track,hw_track,str{tpos}); end

                    % random start parameters in between the given limits
                    % generate charge
                    if Flags.AbsMomentum
                        [param_IR(4:6)]=randt(start,Flags);
                        while param_IR(6) <= pt_limit
                            % this loop is basically obsolete for pt_limit < 0.1 the
                            % while-condtion is then aproximately never fullfilled
                            [param_IR(4:6),start]=randt(start,Flags);
                        end
                    else
                        [param_IR(4:6)]=randt(start,Flags);
                    end

                    % consider charge
                    convf=0.299792458*1e-3*unit;  % conversion factor [Gev/c T^(-1) m^(-1)]
                    q=sign(-0.5+rand(1));
                    convf=convf*q;

                    % Compute curvature from momentum and B field
                    param_IR(6)=convf*Bz/param_IR(6);

                    % Propagation to beamtube
                    [param_start,ierr]=pcyl(radius.b(1),param_IR,sinbmx);
                    switch ierr
                        case 1
                            if octave
                                str=['Error: Vertex outside beamtube!'];
                                fprintf(fidlog,'%s\n',str);
                            else
                                delete(hw_event);delete(hw_track);
                                str=get(whandle,'String');
                                str{2}=['Error: Vertex outside beamtube!'];
                                str{3}=' ';
                                set(whandle,'String',str);
                                set(whandle,'foregroundcolor','red')
                                drawnow
                                fprintf(fidlog,'%s\n',str{2});
                                fprintf(fidlog,'%s\n',str{3});
                            end
                            error('Vertex outside beamtube!');
                            
                        case 2
                            MCpullhit(ITrack)=0;pullhit(ITrack)=0;
                            message=['Track ',num2str(ITrack),...
                                ', propagation to beamtube, max(sin(beta)) exceeded!'];
                            fprintf(fidlog,'%s\n',message);
                    end

                    param_start_all(ITrack,:)=param_start;
                    L=param_IR(4)>[theta1,theta2,theta3,theta4];
                    % to find out the region the track will traverse

                    switch sum(L)
                        case 0  % only forward region, no barrel tracking, use data
                            % of forward region
                            if N.FLayer>1
                                region(ITrack)=1; fnum=fnum+1; Lforward(ITrack)=1;
                                reg='ff'; N.Passive=N.fPassive;
                            else
                                MCpullhit(ITrack)=0;pullhit(ITrack)=0;
                            end
                        case 2  % only barrel region, barrel tracking, use data of
                            % barrel region
                            if N.BLayer>1
                                region(ITrack)=3; bnum=bnum+1; Lbarrel(ITrack)=1;
                                reg='bb'; N.Passive=N.bPassive;
                            else
                                MCpullhit(ITrack)=0;pullhit(ITrack)=0;
                            end
                        case 4  % only rear region, no barrel tracking, use data of
                            % rear region
                            if N.RLayer>1
                                region(ITrack)=1; Lforward(ITrack)=1; rnum=rnum+1;
                                reg='rr'; N.Passive=N.rPassive;
                            else
                                MCpullhit(ITrack)=0;pullhit(ITrack)=0;
                            end
                        otherwise   % track in intermediate region
                            if param_start(3)>=pi/2 reg='ir';
                            else reg='if'; end
                            region(ITrack)=2; inum=inum+1; Linterm(ITrack)=1;
                    end % switch sum(L)

                    %%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if MCpullhit(ITrack)

                        % call of function 'simulation'
                        %disp('simulation')
                        [pstartMS,param_true,coord1,coord2,meas,...
                            hitpattern,reftype,refindex,ref,Xsim,varMSsim,Vsim,...
                            Hsim]=simulation(param_start,region(ITrack));
                        
                        hitpattern_all(ITrack,1:length(hitpattern))=hitpattern;
                        meas_all(ITrack,1:length(meas))=meas;
                        param_true_all(ITrack,:,:)=param_true(1,:);

                        pulldet(ITrack)=min(find(meas>0));
                        if isempty(pulldet(ITrack))
                            pullhit(ITrack)=0;
                        end
                        pulltype(ITrack)=reftype(pulldet(ITrack));
                        pullhit(ITrack)=meas(pulldet(ITrack));

                        L1=(meas==3 | meas==1);
                        L2=(meas==3 | meas==2);

                        if (sum(L1)<minmeas | sum(L2)<minmeas)
                            MCpullhit(ITrack)=0;
                            pullhit(ITrack)=0;
                            if octave
                                str=['Event ',num2str(IEvent),', Track ',...
                                    num2str(ITrack),': Less than ',num2str(minmeas),...
                                    ' measurements, track excluded in computation of pulls!'];
                                fprintf(fidlog,'%s\n',str);
                            else
                                str=get(whandle,'String');
                                str{2}=['Event ',num2str(IEvent),', Track ',...
                                    num2str(ITrack),': Less than ',num2str(minmeas),...
                                    ' measurements'];
                                str{3}=['Track excluded in computation of pulls!'];
                                set(whandle,'String',str);
                                set(whandle,'foregroundcolor','red')
                                drawnow
                                fprintf(fidlog,'%s\n',str{2});
                                fprintf(fidlog,'%s\n',str{3});
                            end
                                
                        end

                        L=meas;
                        L(L==2)=1;
                        L(L==3)=2;
                        ndf(ITrack)=sum(L)-5;   % number of degrees of freedom
                        % of the present track


                        if Flags.Recon

                            % comment this line in to use start parameters without multiple
                            % scattering for the reference track - CAUTION this can cause
                            % errors in the chi^2 tests
                            %pstartMS=param_start;

                            % call of function 'reftrack'
                            %disp('reftrack')
                            [paramr,Hk,Vk,Ak,varMS]=reftrack...
                                (pstartMS,param_start,reftype,refindex,ref);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %disp('reconstruction')
                            [param_fit,Cov,chi2(ITrack)]=reconstruction...
                                (paramr,Hsim,Vk,Ak,coord1,coord2,meas,reftype,...
                                refindex,varMS);
                            
                            param_fit_all(ITrack,:)=param_fit(1,:);
                            Cf_store_all(ITrack,:,:)=Cov(1,:,:);
                            %Cf_store_alldets(ITrack,:,:,:)=Cov;
                            % store fitted parameters at inner side of beamtube and
                            % their corresponding covariance matrix
                            %param_fit_alldets(ITrack,:,:)=param_fit;
                        end % if Flags.Recon

                    end % if MCpullhit(ITrack)

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if MCpullhit(ITrack)
                        if Flags.Recon==1
                            a=[1,pulldet(ITrack)];
                            % to pick only beamtube and first detector layer
                            % call of function 'testquantities'
                            if pulltype(ITrack) Radius=radius.b(a);
                            else                Radius=radius.b(1); end
                            %disp('tests')
                            [mcpull,MC_res(ITrack,:),pull(ITrack,:),res(ITrack,:),...
                                res_true(ITrack,:),MC_chi2(ITrack)]=testquantities...
                                (param_true(a,:),param_fit(a,:),coord1(a(2)),...
                                coord2(a(2)),squeeze(Vk(a(2),:,:)),...
                                squeeze(Hsim(a(2),:,:)),Cov(a,:,:),pullhit(ITrack),...
                                Radius);
                            
                            MC_pull(ITrack,:)=mcpull;
                            L=abs(mcpull)>=limpull;
                            % finds the MC-pulls greater than limpull

                            if ((sum(L)~=0 & Flags.DispBadTrack~=0) | ~isreal(mcpull))
                                % steps in if any pull is greater than limpull
                                % displays pulls and parameters of bad
                                % tracks
                                if octave
                                    fprintf(fidlog,'%s\n',['Bad Track ',num2str(ITrack),' Event ',num2str(IEvent)]);
                                    fprintf(fidlog,'%s\n',['   ','MCpull = ',num2str(mcpull,'% 0.2f')]);
                                    fprintf(fidlog,'%s\n',['   ','Start params: ',num2str(param_start,'% 0.4f')]);
                                else
                                    str=get(whandle,'String');
                                    str{2}=['Bad Track ',num2str(ITrack),' Event ',num2str(IEvent)];
                                    str{3}=['MCpull = ',num2str(mcpull,'% 0.2f')];
                                    str{4}=['Start params: ',num2str(param_start,'% 0.4f')];
                                    set(whandle,'String',str);
                                    set(whandle,'foregroundcolor','red')
                                    drawnow
                                    fprintf(fidlog,'%s\n',str{2});
                                    fprintf(fidlog,'%s\n',['   ',str{3}]);
                                    fprintf(fidlog,'%s\n',['   ',str{4}]);
                                end
                                
                                badcount=badcount+1;
                                % counts the total number of bad tracks
                                bad(L)=bad(L)+1;
                                % increases number of bad pulls for each parameter
                                % seperately
                            end % sum(L)~=0 ...
                        end % if Flags.Recon==1 - end
                    end % if MCpullhit(ITrack) - end
                    if (~octave & interrupt),break,end
                end % end Trackloop

                %--------------------------------------------------------------------------
                %--------------------- Trackloop end --------------------------------------
                %--------------------------------------------------------------------------

                % determines how many tracks of an event are used for the
                % calculation of MC-pulls
                MC_tracknr=sum(MCpullhit);

                if (~octave & interrupt),break,end

                % Collect data from events
                warning off;
                AMCpullhit=logical([AMCpullhit,MCpullhit]);
                AMC_tracknr=[AMC_tracknr,MC_tracknr];
                Apullhit=[Apullhit,pullhit];
                Aparam_start=[Aparam_start;param_start_all];
                ALbarrel=logical([ALbarrel,Lbarrel]);
                Avertex=[Avertex;vertex];
                ALforward=logical([ALforward,Lforward]);
                Aoffset=[Aoffset;offset];
                ALinterm=logical([ALinterm,Linterm]);
                Aparam_true=[Aparam_true;param_true_all];
                Acoord1=[Acoord1;coord1_all];
                Acoord2=[Acoord2;coord2_all];
                Aregion=[Aregion,region];
                Ameas=[Ameas;meas_all];
                %Abadcount=Abadcount+badcount; Abad=Abad+bad;

                if IEvent==1
                    Ahitpattern=[Ahitpattern;hitpattern_all];
                else
                    a=size(hitpattern_all);
                    A=size(Ahitpattern);
                    a=a(2); A=A(2);

                    if a>A
                        Ahitpattern(:,a+1)=Ahitpattern(:,end);
                        Ahitpattern(:,end)=[];
                    elseif a<A
                        hitpattern_all(:,A+1)=hitpattern_all(end);
                        hitpattern_all(:,end)=[];
                    end
                    Ahitpattern=[Ahitpattern;hitpattern_all];
                end

                if Flags.Recon
                    Apulldet=[Apulldet,pulldet];
                    Apulltype=[Apulltype,pulltype];
                    AMC_pull=[AMC_pull;MC_pull];
                    AMC_res=[AMC_res;MC_res];
                    Apull=[Apull;pull];Ares=[Ares;res];%Achi22=[Achi22,chi22];
                    Ares_true=[Ares_true;res_true];Achi2=[Achi2;chi2];
                    ACf_store=[ACf_store;Cf_store_all];Andf=[Andf,ndf];
                    Aparam_fit=[Aparam_fit;param_fit_all];
                    %Aparam_fit_alldets=[Aparam_fit_alldets;param_fit_alldets];
                    %ACf_store_alldets=[ACf_store_alldets;Cf_store_alldets];
                    Aresout_coordinate=[Aresout_coordinate,resout_coordinate];
                    %Achi_hist=[Achi_hist,chi_hist];
                    AMC_chi2=[AMC_chi2,MC_chi2];
                end

            end % end Eventloop

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%% Eventloop end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ~octave & interrupt
                delete(hw_track)
                if Flags.Curves
                    if N.Event==1
                        message=['Interrupted in curve ',num2str(ICurve),', point ',num2str(IPoint),...
                            ', track ',num2str(ITrack),'!'];
                    else
                        message=['Interrupted in curve ',num2str(ICurve),', point ',num2str(IPoint),...
                            ', event ',num2str(IEvent),', track ',num2str(ITrack),'!'];
                        delete(hw_event);
                    end
                else
                    if N.Event==1
                        message=['Interrupted in track ',num2str(ITrack),'!'];
                    else
                        message=['Interrupted in event ',num2str(IEvent),' track ',num2str(ITrack),'!'];
                        delete(hw_event);
                    end
                end
                fprintf(fidlog,'%s\n',message); fclose(fidlog);
                str=get(whandle,'String');
                str{2}=message;
                str{3}=[];
                str{4}=[];
                set(whandle,'String',str);
                set(whandle,'foregroundcolor','red');
                drawnow
                
                %warndlg(str,'Simulation interrupted');
                return
            end

            if Flags.wait
                if N.Event>1 hwait=hw_event; delete(hw_track);
                else         hwait=hw_track; end
            end

            if Flags.Recon
                [rms,imp]=logfile(start,Aresout_coordinate,Aoffset,Apull,Apullhit,AMCpullhit,...
                    AMC_pull,Apulltype,bnum,fnum,rnum,inum,bad,...
                    badcount,AMC_res,Ares_true,Aparam_start,AMC_tracknr,...
                    Avertex,Aparam_fit,Achi2,Andf,AMC_chi2,ACf_store,Aregion,...
                    barrelfile,forwardfile);
                %----------------------for calculation of curves------
                %rms(1)=RPhi rms(2)=z rms(3)=theta rms(4)=phi rms(5)=dpt/pt rms(6)=dpt/pt^2
                %imp(1)=mean(ip2d), imp(2)=std(ip2d), imp(3)=mean(ip3d), imp(4)=std(ip3d)
                %
                resultsofcurves(1,ICurve,IPoint)=rms(1);
                resultsofcurves(2,ICurve,IPoint)=rms(2);
                resultsofcurves(3,ICurve,IPoint)=rms(3);
                resultsofcurves(4,ICurve,IPoint)=rms(4);
                resultsofcurves(5,ICurve,IPoint)=rms(5);
                resultsofcurves(6,ICurve,IPoint)=rms(6);
                resultsofcurves(7,ICurve,IPoint)=imp(2);
                resultsofcurves(8,ICurve,IPoint)=imp(3);
                %}---------------------------------------------

            end % Flags.Recon
        end % if Flags.Simul==1 - end

        %------------ Close log file -----------------------------------------------
        %fclose(fidlog);
        %
        %R=radius.b(1);
        %a=AMCpullhit;
        %param_start=Aparam_start(a,:); param_fit=Aparam_fit(a,:);
        %Cf_store=ACf_store(a,:,:); pulltype=Apulltype(a);
        %ndf=Andf(a); chi2=Achi2(a); pullhit=Apullhit(a);
        %resout_coordinate=Aresout_coordinate(a);
        %save results AMC_pull Apull Ahitpattern ACf_store  ...
        %    Aparam_start Ameas vertex param_start param_fit R Cf_store ...
        %    Bz unit Flags fidlog hwait AMC_tracknr Avertex SPR...
        %    Cf_store Aoffset pulltype ndf chi2 pullhit resout_coordinate...
        %    convf N GeomVersion runnumber Flags Apullhit Apulltype AMCpullhit...
        %    AMC_res Ares_true Aparam_fit limpull;
        %save rave vertex param_start param_fit R Cf_store ...
        %    Bz unit Flags fidlog hwait;
        %save JAS3 AMC_tracknr R Avertex param_start param_fit SPR...
        %    Cf_store Aoffset pulltype ndf chi2 pullhit resout_coordinate...
        %    Bz unit convf N GeomVersion runnumber Flags
        %save histograms AMC_pull Apull Apullhit Apulltype AMCpullhit...
        %    SPR R AMC_res Ares_true Aparam_start Aparam_fit convf unit Bz limpull
        %}
        if Flags.Simul & Flags.wait delete(hwait); end

        %---------------for calculation of curves-------------------------
        %
    end % end pointloop
end % end curveloop
R=radius.b(1);
a=AMCpullhit;
param_start=Aparam_start(a,:); param_fit=Aparam_fit(a,:);
Cf_store=ACf_store(a,:,:); pulltype=Apulltype(a);
ndf=Andf(a); chi2=Achi2(a); pullhit=Apullhit(a);
resout_coordinate=Aresout_coordinate(a);
if Flags.Curves
    if octave
        save('-v7','results.mat');
    else
        save('results.mat');
    end
    %save('results.mat','resultsofcurves','AMC_pull','Apull','Ahitpattern','ACf_store',...
    %        'Aparam_start','Ameas','vertex','param_start','param_fit','R','Cf_store',...
    %        'Bz','unit','Flags','fidlog','hwait','AMC_tracknr','Avertex','SPR'...
    %        'Cf_store','Aoffset','pulltype','ndf','chi2','pullhit','resout_coordinate',...
    %        'convf','N','GeomVersion','runnumber','Flags','Apullhit','Apulltype','AMCpullhit',...
    %        'AMC_res','Ares_true','Aparam_fit','limpull','nameb','namef');
    LDT_curves('results.mat');
    %
    %quantity={'rmsRPhi','rmsz','rmstheta','rmsphi','rmsdptpt','rmsdptpt2','stdimp2d','meanimp3d'};
    %use=[6,7];
    %for d=use
        %str=[forwardfile(6:end),'-',barrelfile(6:end),'-',quantity{d},'.mat'];
    %    str=['sine_LDCprime/LDCprime-zover100-',quantity{d},'.mat'];
    %    fid=fopen(str,'w');
    %    formatstr=[];
    %    for k=1:N.Point
    %        formatstr=[formatstr,'%g   '];
    %    end
    %    formatstr=[formatstr,'\n'];
    %    %fprintf(fidlog,formatstr,dptpt2');
    %    fprintf(fid,formatstr,squeeze(resultsofcurves(d,:,:))');
    %    fclose(fid);
    %end
    %}
else
    if octave
        save('-v7','results.mat');
    else
        save('results.mat');
    end
    %save('results.mat');
    %save results AMC_pull Apull Ahitpattern ACf_store  ...
    %    Aparam_start Ameas vertex param_start param_fit R Cf_store ...
    %    Bz unit Flags fidlog hwait AMC_tracknr Avertex SPR...
    %    Cf_store Aoffset pulltype ndf chi2 pullhit resout_coordinate...
    %    convf N GeomVersion runnumber Flags Apullhit Apulltype AMCpullhit...
    %    AMC_res Ares_true Aparam_fit limpull;
    LDT_residuals('results.mat');
end
%toc
t=toc;
strfinished=['Simulation successfully finished in ',num2str(t),' seconds!'];
fprintf(fidlog,'\n%s\n',strfinished);
fclose(fidlog);

if octave
    %str=['\nSimulation successfully finished\n\n'];
    printf(strfinished);
    fflush(stdout);
else
    str=get(mhandle,'String');
    %str{2}=['Simulation successfully finished in ',num2str(t),' seconds!'];
    str{2}=strfinished;
    str{3}=[];
    str{4}=[];
    set(mhandle,'String',str);
    drawnow
end
   %}
%----------------------------------------------------------------