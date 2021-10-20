function LDT_residuals(fnameresults)

% function HISTOGRAMS
% Called by output
% Main program: LDT_main
%
% Input:    MC_pull     Array of pull quantities, computed from the start
%                       parameters and the fitted parameters at the inner
%                       side of the beamtube, for every track and every
%                       coordinate (Phi,z,theta,beta,kappa)
%                       (MC stands for Monte Carlo)
%           MCpull_mean Mean of MC_pull over all tracks, for every coord.
%           MCpull_std  Standard deviation of MC_pull over all tracks, for
%                       every coordinate
%           pull        Real pull quantities at the first detector layer,
%                       computed from the fitted 'measurements' and the
%                       real (simulated) measurements, for every track and
%                       coordinate ( (RPhi,z') and (u,v) respectively)
%           bpull_mean  Mean of pull over all tracks from barrel region
%           bpull_std   Standard deviation of pull over all barrel tracks
%           fpull_mean  Mean of pull over all tracks from forward region
%           fpull_std   Standard deviation of pull over all forward tracks
%	    pulltype	Array, which indicates if pulls for the forward or
%			barrel region should be histogramed
%           pullhit     Logical array, which indicates the tracks for those
%                       real pulls can be computed
%           MCpullhit   Logical array, which indicates the tracks for those
%                       Monte-Carlo pulls can be computed
%           MC_hist     Arrays containing the data for MC-residuals histograms
%
% Output: none
%
%   HISTOGRAMS displays the pull quantities at the inner side of the
%   beamtube and at the innermost detector together in one figure.
%   Furthermore it creates one more figure and draws histograms of the
%   MC-residuals of RPhi, z, theta, phi, dpt/pt and dpt/pt^2
%

%global convf unit
%global Bz
%global limpull
%global Flags
warning off;
close(figure(3))
%action='residuals';

if nargin==0
    fnameresults='results.mat';
end
%(MC_pull,MCpull_mean,MCpull_std,pull,pullhit,pulltype,MCpullhit,SPR,rms,MC_hist)
load(fnameresults,'AMC_pull','Apull','Apullhit','Apulltype','AMCpullhit',...
            'SPR','R','AMC_res','Ares_true','Aparam_start','Aparam_fit','convf','unit','Bz','limpull');
MC_pull=AMC_pull;pull=Apull;pullhit=Apullhit;pulltype=Apulltype;MCpullhit=AMCpullhit;
%[rms,MC_hist]=mcrms_test(R,AMC_res,Ares_true,Aparam_start,Aparam_fit,AMCpullhit);
[rms,MC_hist]=mcrms(R,AMC_res,Ares_true,Aparam_start,Aparam_fit,AMCpullhit);
Luse=(AMCpullhit>0); % picks useful tracks
if sum(Luse)~=0
    MCpull_mean=squeeze(mean(AMC_pull(AMCpullhit,:))); % mean of MC-pulls
    MCpull_std=squeeze(std(AMC_pull(AMCpullhit,:)));   % std dev. of MC-pulls
else
    MCpull_mean=nan*ones(1,5);
    MCpull_std=nan*ones(1,5);
end
format short e
dots=200;   % number of dots used to draw the lines


%------------------------Results - Residuals-------------------------------

% Histograms of results
L=MCpullhit;
numtrack=sum(L);

figure(3) % opens a figure
figure(3)
set(gcf,'Position',[95,50,900,650])
set(gcf,'Name',['Residuals (',fnameresults,')']);
subplot(2,4,6)
% pick region in the midth of the second row for
% explaining text
hold on
axis off
%h=text(0,0.7,{'Top:','',''});
h=text(0,0.7,'Top:');
set(h,'fontweight','bold');
%text(0.4,0.7,{'Residuals at the beamtube,','difference of true and fitted','track parameters'});
text(0.5,0.7,'Residuals at the beamtube,');
text(0.5,0.6,'difference of true and fitted');
text(0.5,0.5,'track parameters');
h=text(0,0.4,'Left:');
set(h,'fontweight','bold');
text(0.5,0.38,'Relative deviation \Delta p_t/p_t');
h=text(0,0.2,'Right:');
set(h,'fontweight','bold');
text(0.5,0.19,'Relative deviation \Delta p_t/p_t^2');
%hold off
% position, width, length
%subplot(2,4,6)
% pick region in the midth of the second row for
% explaining text
%hold on
%axis off
text(0 ,1  ,[num2str(numtrack),' tracks with ']);
text(1 ,1  ,[num2str(SPR.Thetamin(end)),' <= \vartheta <= ',num2str(SPR.Thetamax(end)),' deg']);
text(1 ,0.9,[num2str(SPR.Ptmin(end)),' <= p_t <= ',num2str(SPR.Ptmax(end)),' GeV']);
hold off

% Residuals of RPhi at 1st detector
k=0;
subplot(2,4,1)                     % pick upper left corner

lim=max(abs(MC_hist(:,1)));
% Protect against no useful tracks
if isempty(lim)
    lim=1;
end
% finds out maximal value to fill into histogram
hist(MC_hist(:,1),linspace(-lim,lim));  % draws histogram
n=hist(MC_hist(:,1),linspace(-lim,lim));% n holds bin contents
xl=max(abs(get(gca,'xlim')));
xlim([-xl xl]);

%title(['r.m.s.=',num2str(rms(1)*1000*unit,'%6.4g'),' \mu{}m']);
xlabel(['r.m.s.=',num2str(rms(1)*1000*unit,'%6.4g'),' \mu{}m']);
ylabel('\Delta(R\Phi)');
hold on
n=max(n);
% maximal content of a bin, used to draw the lines of rms
plot(-rms(1),linspace(0,n,200),'r');
% plots root mean square as a red line
plot(rms(1),linspace(0,n,200),'r');
hold off

% Residuals of phi (=beta+Phi) at 1st detector
subplot(2,4,2)  % pick upper right corner

lim=max(abs(MC_hist(:,4)));
% Protect against no useful tracks
if isempty(lim)
    lim=1;
end
hist(MC_hist(:,4),linspace(-lim,lim));
xl=max(abs(get(gca,'xlim')));
xlim([-xl xl]);
ylims=get(gca,'ylim');
xlabel(['r.m.s.=',num2str(rms(4)*1000,'%6.4g'),' mrad']);
ylabel('\Delta(\phi)');
hold on
plot([-rms(4) -rms(4)],ylims,'r');
plot([rms(4) rms(4)],ylims,'r');
hold off


for k=2:3   % draw histograms of z and theta in a loop
    % Residuals of z and theta at 1st detector
    subplot(2,4,k+1)  % pick 2nd and 3rd position of upper row
    lim=max(abs(MC_hist(:,k)));
    % Protect against no useful tracks
    if isempty(lim)
        lim=1;
    end
    hist(MC_hist(:,k),linspace(-lim,lim));
    xl=max(abs(get(gca,'xlim')));
    xlim([-xl xl]);
    ylims=get(gca,'ylim');
    %xlabel('RMS=',num2str(rms(k))});
    hold on
    plot([-rms(k) -rms(k)],ylims,'r');
    plot([rms(k) rms(k)],ylims,'r');
    hold off

    switch k    % different titles
        case 2
            %title({'MC-Residuals','of z'});
            ylabel('\Delta(z)');
            xlabel(['r.m.s.=',num2str(rms(k)*1000*unit,'%6.4g'),' \mu{}m']);
        otherwise
            %title({'MC-Residuals','of \vartheta'});
            ylabel('\Delta(\vartheta)');
            xlabel(['r.m.s.=',num2str(rms(k)*1000,'%6.4g'),' mrad']);
    end
end % for k=2:3 - end



% relative deviation delta pt/pt
subplot(2,4,5)  % pick lower left corner
lim=max(abs(MC_hist(:,5)));
% Protect against no useful tracks
if isempty(lim)
    lim=1;
end
L(abs(MC_hist(:,5))>lim)=0;
hist(MC_hist(:,5),linspace(-lim,lim));
ylabel('\Delta p_t/p_t');
xl=max(abs(get(gca,'xlim')));
xlim([-xl xl]);
ylims=get(gca,'ylim');
xlabel(['r.m.s.=',num2str(rms(5),'%6.4E')]);
hold on
plot([-rms(5) -rms(5)],ylims,'r');
plot([rms(5) rms(5)],ylims,'r');
hold off

% relative deviation delta pt/pt^2
subplot(2,4,8)  % pick lower right corner
lim=max(abs(MC_hist(:,6)));
% Protect against no useful tracks
if isempty(lim)
    lim=1;
end
hist(MC_hist(:,6),linspace(-lim,lim));
ylabel('\Delta p_t/p_t^2');
xl=max(abs(get(gca,'xlim')));
xlim([-xl xl]);
ylims=get(gca,'ylim');
xlabel(['r.m.s.=',num2str(rms(6),'%6.4E')]);
hold on
plot([-rms(6) -rms(6)],ylims,'r');
plot([rms(6) rms(6)],ylims,'r');
hold off