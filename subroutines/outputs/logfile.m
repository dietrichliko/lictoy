function [rms,imp]=logfile(SPR,Aresout_coordinate,Aoffset,Apull,Apullhit,AMCpullhit,AMC_pull,Apulltype,bnum,fnum,rnum,inum,bad,...
            badcount,AMC_res,Ares_true,Aparam_start_all,AMC_tracknr,...
            Avertex,Aparam_fit_all,Achi2,Andf,AMC_chi2,ACf_store_all,Aregion,barrelfile,forwardfile)
            
% function logfile
% Called by: LDT_main
% Main function: LDT_main
%
% LOGFILE writes the simulation parameters, the test results and the reconstruction 
% results to the log file.
        
global Flags Bz Mass N fidlog limpull unit radius hwait mhandle
warning off;

Luse=(AMCpullhit>0); % picks useful tracks
if sum(Luse)~=0
    MCpull_mean=squeeze(mean(AMC_pull(AMCpullhit,:))); % mean of MC-pulls
    MCpull_std=squeeze(std(AMC_pull(AMCpullhit,:)));   % std dev. of MC-pulls
else
    MCpull_mean=nan*ones(1,5);
    MCpull_std=nan*ones(1,5);
end
Lb=( (Apullhit>0) & Aregion==3 ) ; % picks useful tracks from barrel
if sum(Lb)~=0
    bpull_mean(1)=mean(Apull(Lb&Apullhit~=2,1)); % mean of pulls (1st coord)
    bpull_std(1)=std(Apull(Lb&Apullhit~=2,1));   % standard deviation of pulls (1st coord)
    bpull_mean(2)=mean(Apull(Lb&Apullhit~=1,2)); % mean of pulls (2nd coord)
    bpull_std(2)=std(Apull(Lb&Apullhit~=1,2));   % standard deviation of pulls (2nd coord)
else
    %bpull_mean=[NaN,NaN];
    %bpull_std=[NaN,NaN];
    bpull_mean=[0,0];
    bpull_std=[0,0];
end
Lf=( (Apullhit>0) & Aregion==1 ); % picks useful tracks from forward
if sum(Lf)~=0
    fpull_mean(1)=mean(Apull(Lf&Apullhit~=2,1)); % mean of pulls (1st coord)
    fpull_std(1)=std(Apull(Lf&Apullhit~=2,1));% std dev of pulls (1st coord)
    fpull_mean(2)=mean(Apull(Lf&Apullhit~=1,2)); % mean of pulls (2nd coord)
    fpull_std(2)=std(Apull(Lf&Apullhit~=1,2));% std dev of pulls (2nd coord)
else
    fpull_mean=[0,0];
    fpull_std=[0,0];
end
%
Li=( (Apullhit>0) & Aregion==2 & Apulltype==1); % picks useful tracks from intermediate (barrel pull)
if sum(Li)~=0
    ibpull_mean(1)=mean(Apull(Li&Apullhit~=2,1)); % mean of pulls (1st coord)
    ibpull_std(1)=std(Apull(Li&Apullhit~=2,1));% std dev of pulls (1st coord)
    ibpull_mean(2)=mean(Apull(Li&Apullhit~=1,2)); % mean of pulls (2nd coord)
    ibpull_std(2)=std(Apull(Li&Apullhit~=1,2));% std dev of pulls (2nd coord)
else
    ibpull_mean=[0,0];
    ibpull_std=[0,0];
end

Li=( (Apullhit>0) & Aregion==2 & Apulltype==0); % picks useful tracks from intermediate (forward pull)
if sum(Li)~=0
    ifpull_mean(1)=mean(Apull(Li&Apullhit~=2,1)); % mean of pulls (1st coord)
    ifpull_std(1)=std(Apull(Li&Apullhit~=2,1));% std dev of pulls (1st coord)
    ifpull_mean(2)=mean(Apull(Li&Apullhit~=1,2)); % mean of pulls (2nd coord)
    ifpull_std(2)=std(Apull(Li&Apullhit~=1,2));% std dev of pulls (2nd coord)
else
    ifpull_mean=[0,0];
    ifpull_std=[0,0];
end
%}

%------------------- Output of simulation parameters and results ------------------------

% title line
fprintf(fidlog,'%s\n','-----------------------------------------------------------------');
fprintf(fidlog,'%s\n','-------------------   SIMULATION PARAMETERS   -------------------');
fprintf(fidlog,'%s\n','-----------------------------------------------------------------');
%if Flags.AbsMomentum
%    fprintf(fidlog,'%s %8.4f %8.4f %s\n','P_abs:[Gev] [',SPR.P_abs,']');
%else
    fprintf(fidlog,'%s %8.4f %8.4f %s\n','Pt:[Gev]    [',SPR.Pt,']');
%end
fprintf(fidlog,'%s %8.4f %8.4f %s\n','phi:[rad]   [',SPR.phi,']');
fprintf(fidlog,'%s %8.4f %8.4f %s\n','Theta:[rad] [',SPR.Theta,']');
fprintf(fidlog,'%s %8.4f %8.4f %s\n','x:[mm]      [',SPR.x,']');
fprintf(fidlog,'%s %8.4f %8.4f %s\n','y:[mm]      [',SPR.y,']');
fprintf(fidlog,'%s %8.4f %8.4f %s\n\n','z:[mm]      [',SPR.z,']');
fprintf(fidlog,'%s\t %8.4f %s\n','Magnetic Field Bz:   ',Bz,' [T]');
fprintf(fidlog,'%s\t %8.4f %s\n','Mass of the particles',Mass,' [GeV]');
fprintf(fidlog,'\t%s\n','(plays only a role in multiple scattering =>');
fprintf(fidlog,'\t%s\n','zero mass means all particles are treated relativisticly!)');

fprintf(fidlog,'\n%s\n','Flags:');
if Flags.ThetaSymmetry
    fprintf(fidlog,'%s\n','Symmetry in Theta         : on');
    fprintf(fidlog,'%s\n','Theta range in two slices (forward and rear)');
else
    fprintf(fidlog,'%s\n','Symmetry in Theta         : off');
end
if Flags.AbsMomentum
    fprintf(fidlog,'%s\n','Use absolute momentum     : on');
else
    fprintf(fidlog,'%s\n','Use absolute momentum     : off');
end
if Flags.ScaleDownTPC
    fprintf(fidlog,'%s\n','Scale down TPC by factor 5: on');
else
    fprintf(fidlog,'%s\n','Scale down TPC by factor 5: off');
end
if Flags.MulSca
    fprintf(fidlog,'%s\n','Multiple scattering       : on');
else
    fprintf(fidlog,'%s\n','Multiple scattering       : off');
end
if Flags.Smear
    fprintf(fidlog,'%s\n','Measurement errors        : on');
else
    fprintf(fidlog,'%s\n','Measurement errors        : off');
end
if Flags.DispBadTrack
    fprintf(fidlog,'%s\n','Display bad tracks        : on');
else
    fprintf(fidlog,'%s\n','Display bad tracks        : off');
end
if Flags.Chi2
    fprintf(fidlog,'%s\n','Chi2                      : on');
else
    fprintf(fidlog,'%s\n','Chi2                      : off');
end
    
fprintf(fidlog,'\n%s\n',['Barrel geometry:  ',barrelfile]);
fprintf(fidlog,'%s\n',['Forward geometry: ',forwardfile]);
fprintf(fidlog,'%s\n\n','-----------------------------------------------------------------');

fprintf(fidlog,'%s\t\t %i\n',  'Number of events:      ',N.Event);
fprintf(fidlog,'%s\t\t %i\n','Number of tracks/event:',N.Track);

fprintf(fidlog,'\t%s\t %i\n','barrel region:        ',bnum);
fprintf(fidlog,'\t%s\t %i\n','forward region:       ',fnum);
fprintf(fidlog,'\t%s\t %i\n','rear region:          ',rnum);
fprintf(fidlog,'\t%s\t %i\n\n','intermediate region:  ',inum);
fprintf(fidlog,'%s\n','-----------------------------------------------------------------');
fprintf(fidlog,'%s\n','-------------------------    RESULTS   --------------------------');
fprintf(fidlog,'%s\n','-----------------------------------------------------------------');

str=[num2str(sum(AMCpullhit)),' out of ',...
    num2str(length(AMCpullhit)),' tracks are used for test statistics'];
fprintf(fidlog,'%s\n\n',str);

if Flags.DispBadTrack==1
    % displays the value of limpull
    str=['Pulls greater than ',num2str(limpull),...
        ' are considered bad'];
    fprintf(fidlog,'%s\n\n',str);
    % displays the total number of bad pulls for every coordinate
    % (Phi,z,theta,beta,1/R)
    str='Bad pulls:';
    fprintf(fidlog,'%s\n',str);
    str=num2str(bad);
    fprintf(fidlog,'%s\n\n',str);
    % displays the total number of bad tracks
    str='Total number of bad tracks:';
    fprintf(fidlog,'%s\n',str);
    str=num2str(badcount);
    fprintf(fidlog,'%s\n\n',str);
end
% displays the means of Monte-Carlo pulls for every coordinate
% (Phi,z,theta,beta,1/R)
format +
format short
%dig=4;
str=['Standardized residuals'];
fprintf(fidlog,'%s\n',str);
str='      Phi            z          theta        beta       kappa';
fprintf(fidlog,'%s\n',str);
str=['mean: ',num2str(MCpull_mean)];
fprintf(fidlog,'%s\n',str);
str=['std:  ',num2str(MCpull_std)];
fprintf(fidlog,'%s\n\n',str);
% displays the means of pulls at innermost detector layer,
% for every coordinate (RPhi,z',u,v)
str=['Pulls of innermost hit'];
fprintf(fidlog,'%s\n',str);
str=('      Barrel tracks  Fwd/Rear tracks  -----Intermediate tracks-----');
fprintf(fidlog,'%s\n',str);
str=('       RPhi    z       u       v       RPhi    z       u       v');
fprintf(fidlog,'%s\n',str);
str=(['mean: ',sprintf('%+0.2f   ',bpull_mean),sprintf('%+0.2f   ',fpull_mean),sprintf('%+0.2f   ',ibpull_mean),sprintf('%+0.2f   ',ifpull_mean)]);
fprintf(fidlog,'%s\n',str);
str=(['std:  ',sprintf('%+0.2f   ',bpull_std),sprintf('%+0.2f   ',fpull_std),sprintf('%+0.2f   ',ibpull_std),sprintf('%+0.2f   ',ifpull_std)]);
fprintf(fidlog,'%s\n\n',str);

% output of the root mean square of the MC-residuals into the log
% file

[rms,mc_hist]=mcrms(radius.b(1),AMC_res,Ares_true,Aparam_start_all,Aparam_fit_all,AMCpullhit);
%[rms,mc_hist]=mcrms_test(radius.b(1),AMC_res,Ares_true,Aparam_start_all,Aparam_fit_all,AMCpullhit);
fprintf(fidlog,'%s \n','R.M.S. values of the residuals and momenta [mu] and [mrad] resp.');
fprintf(fidlog,'%s \n','      RPhi          z      theta        phi     dpt/pt   dpt/pt^2');
%conversion to [mu*rad,mu,mrad,mrad,nounit,c/GeV]
rms1=[rms(1:2)*1000*unit,rms(3:4)*1000,rms(5:6)];
[str,errm]=sprintf('%10.5f ',rms1);
fprintf(fidlog,'%s \n\n',str);


% calculation of the impact parameters (projected and 3d)
% to hand over only data for which MC-pulls were  calculated
a=AMCpullhit;
Cf=ACf_store_all(a,:,:);
[ip2d,ipm2d,ipstd2d,ip3d,ipm3d,ipstd3d,ipz,ipmz,ipsdtz,deltazaxis]=...
    impact_par(AMC_tracknr,radius.b(1),Avertex,Aparam_fit_all(a,:),Cf,N);
imp=[ipm2d,ipstd2d,ipm3d,ipstd3d];

%
%write the impact parameters into the log file
fprintf(fidlog,'%s \n','Impact paramters [mu]');
fprintf(fidlog,'%s \n','          2d    3d (abs)');
fprintf(fidlog,'%s %8.4f %8.4f\n','mean:',ipm2d*1000*unit,ipm3d*1000*unit);
fprintf(fidlog,'%s %8.4f %8.4f\n\n','std: ',ipstd2d*1000*unit,ipstd3d*1000*unit);
%}

if Flags.Chi2==1

    if sum(Luse)~=0
        chi2_mean=mean(Achi2(AMCpullhit));
        chi2_var=var(Achi2(AMCpullhit));
        ndfmean=sum(Andf(AMCpullhit))/sum(AMCpullhit);
    else
        chi2_mean=nan;
        chi2_var=nan;
        ndfmean=nan;
    end
    %bndf=(N.BLayer-N.bPassive)*2-5;
    % increases number of bad pulls for each parameter seperately
    %fndf=((N.FLayer+N.RLayer)/2-(N.rPassive+N.rPassive)/2)*2-5;
    % displays the mean of the total chi2s
    str=('Chi^2');
    fprintf(fidlog,'%s\n',str);
    % displays number of degrees of freedom
    str=(['ndf:    ',num2str(ndfmean),', average']);
    fprintf(fidlog,'%s\n',str);
    str=(['mean:   ',num2str(chi2_mean)]);
    fprintf(fidlog,'%s\n',str);
    fprintf(fidlog,'%s %f\n\n','ratio: ',ndfmean/chi2_mean);

    fprintf(fidlog,'%s\n','Mean of the MC-Chi^2 at the beamtube:');
    fprintf(fidlog,'%f\n\n',mean(AMC_chi2));
end


%if (Flags.DispPull==1 | Flags.DispRes==1 | Flags.DispImp==1)
    % call of function 'histograms'
%    if Flags.wait
%        waitbar(1,hwait,'Generating histograms');
%    end
    %histograms(AMC_pull,AMC_res,MCpull_mean,MCpull_std,Apull,Ares,...
    %    bpull_mean,bpull_std,fpull_mean,fpull_std,Ares_true,...
    %    Aparam_start_all,Radius(1),Apullhit,AMCpullhit,ALbarrel,...
    %    ALforward,ALinterm,SPR);
    %histograms_test(AMC_pull,MCpull_mean,MCpull_std,Apull,...
    %            Apullhit,Apulltype,AMCpullhit,SPR,rms,mc_hist);

    
    %if sum(Achi22)/length(Achi22) ~= [-1 -1]
    %    histogram(Achi_hist,'firt layer in the TPC','last hit plane detector',5);
    %    histogram(Achi22,'before','after',6);
    %end
%end
