%function resout(impact2d,impact3d,SPR,tracknr,radius,vertex,param_start,param_fit,Cf_store,offset...
%    ,pulltype,ndf,chi2,pullhit,coordinate)
function LDT2jas3(filename,varargin)

% function jas3
% Called by LDT_main
% Main program: LDT_main
%
% JAS3 writes the start parameters, the fitted pararmeters, the 
% covariance matrix at the inner side of the beamtube and the vertex 
% coordinates to the JAS3 input file. The filename is specified in
% momentumfilename. If momentumfilename is empty, a filename is 
% generated using the current date and time


%global Bz unit convf
%global N GeomVersion runnumber
%global Flags mhandle momentumfilename

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
load(fnameresults,'AMC_tracknr','R','Avertex','param_start','param_fit','SPR',...
            'Cf_store','Aoffset','pulltype','ndf','chi2','pullhit','resout_coordinate',...
            'Bz','unit','convf','N','GeomVersion','runnumber','Flags');
tracknr=AMC_tracknr; radius=R; vertex=Avertex; offset=Aoffset; coordinate=resout_coordinate;
[impact2d,ipm2d,ipstd2d,impact3d,ipm3d,ipstd3d,ipz,ipmz,ipsdtz,deltazaxis]=...
    impact_par(tracknr,radius,vertex,param_fit,Cf_store,N);

%if ~isempty(momentumfilename)
%    output=momentumfilename;
if isempty(filename)
    T=clock;
    filename=['track_data_',num2str(T(3)),'.',num2str(T(2)),'.',num2str(T(1)),'_',num2str(T(4)),':',num2str(T(5)),':',num2str(floor(T(6))),'.txt'];
end
seperator=...
    '-----------------------------------------------------------------';

fout=fopen(filename,'w');
disp(['Writing file: ',filename]);
%----------------------------- writing  header ----------------------------

radstr='[rad]';
gevstr='[GeV]';
if unit==1
    ustr='[mm]';
else
    ustr='[m]';
end

% appending 2 zeros for some reason to the version numbers of the input
% sheets
fprintf(fout,'%s\t %s\t %06i\t %i\n','Version: ',GeomVersion.Bfilename,GeomVersion.BNo,runnumber);
fprintf(fout,'%s\t %s\t %06i\t %i\n','Version: ',GeomVersion.Ffilename,GeomVersion.FNo,runnumber);

fprintf(fout,'%s \t%i\n','Number of events:             ',N.Event);
fprintf(fout,'%s \t%i\n','Number of track/event:        ',N.Track);
fprintf(fout,'%s \t%f\n',['Radius of the beam tube ',ustr,': '],radius);
fprintf(fout,'%s\n','Simulation parameters:');

if Flags.AbsMomentum
    fprintf(fout,'%s %8.4f %8.4f %s\n',['P_abs [rad]:     ['],num2str([SPR.P_absmin,SPR.P_absmax]),']');
else
    fprintf(fout,'%s %8.4f %8.4f %s\n',['Pt [GeV]:        ['],num2str([SPR.Ptmin,SPR.Ptmax]),']');
end

fprintf(fout,'%s %8.4f %8.4f %s\n',['phi [rad]:       ['],[SPR.phimin,SPR.phimax],']');
fprintf(fout,'%s %8.4f %8.4f %s\n',['Theta [rad]:     ['],[SPR.Thetamin,SPR.Thetamax],']');
fprintf(fout,'%s %8.4f %8.4f %s\n',['x ',ustr,':          ['],SPR.x,']');
fprintf(fout,'%s %8.4f %8.4f %s\n',['y ',ustr,':          ['],SPR.y,']');
fprintf(fout,'%s %8.4f %8.4f %s\n\n',['z ',ustr,':          ['],SPR.z,']');
fprintf(fout,'%s\t %8.4f\n','Magnetic Field Bz [T]:',Bz);

if ( Flags.ThetaSymmetry )
    fprintf(fout,'%s\n','Theta range in two slices (forward and rear)');
end
fprintf(fout,'%s\n',seperator);
fprintf(fout,'%s\n','Numbers of tracks per vertex:');
[str,errm]=sprintf('%i ',tracknr);
fprintf(fout,'%s\n',str);
fprintf(fout,'%s\n',seperator);
fprintf(fout,'%s\n',['units are ',ustr,' and ',radstr]);

%---------------------------- end of header -------------------------------

index=1;
for eventnr=1:N.Event
    [str,errm]=sprintf('%8.4e ',vertex(eventnr,:));
    fprintf(fout,'\n%s\t\t\t\t %s\n','vertex:',str);
    for k0=1:tracknr(eventnr)
                
        [str,errm]=sprintf('%8.4e ',[param_start(index,:),abs((convf*Bz)/param_start(index,5))]);
        fprintf(fout,'%s\t %s\n','true parameters:',str);
        [str,errm]=sprintf('%8.4e ',param_fit(index,:));            
        fprintf(fout,'%s\t %s\n','fitted parameters:',str);
        [str,errm]=sprintf('%8.4e ',impact2d(index));            
        fprintf(fout,'%s\t %s\n','impact 2d: ',str);
        [str,errm]=sprintf('%8.4e ',impact3d(index));            
        fprintf(fout,'%s\t %s\n','impact 3d: ',str);
        
        
        % writing ndf,chi^2, detector type, coordinate and number of
        % measured coordinates
        fprintf(fout,'%s\t\t\t','ndf - chi2:');
        [str,errm]=sprintf('%6.2f %6.2f',[ndf(index),chi2(index)]);
        fprintf(fout,'%s\n',str);
        fprintf(fout,'%s\t\t %s\t','measurements:',char(pulltype(index)));
        [str,errm]=sprintf('%8.4f %6.2f\n',[coordinate(index),pullhit(index)]);
        fprintf(fout,'%s',str); 
        
        % writing cov-matrix
        triu_cov=triu(squeeze(Cf_store(index,:,:)));
        vec_cov=[triu_cov(1,:),triu_cov(2,2:5),triu_cov(3,3:5),triu_cov(4,4:5),triu_cov(5,5),];
        [str,errm]=sprintf('%8.4e ',vec_cov);
        fprintf(fout,'%s\t %s\n','covariance matrix:',str);
        index=index+1;
       
    end
end

fprintf(fout,'\n%s\n','EOF');
fclose(fout);
disp('finished!');
