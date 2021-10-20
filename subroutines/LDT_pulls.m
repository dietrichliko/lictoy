function LDT_pulls(fnameresults)

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
%action='pulls';

warning off;
close(figure(2));

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

%----------------------Test Pulls------------------------------------------


% Histograms of test quantities
add=1;

h=figure(2);   % test quantities are displayed in figure 2
figure(2)
set(h,'Position',[60,65,900,650])    % position, width, length
set(h,'Name',['Pull quantities (',fnameresults,')']);

% RPhi-pulls at 1st detector
subplot(2,5,1)  % pick upper left corner
box on;
L=(pullhit>0)&pulltype==1;
if sum(L)~=0
    % means and standard deviations
    bpull_mean(1)=mean(pull(L&pullhit~=2,1)); % mean of pulls (1st coord)
    bpull_std(1)=std(pull(L&pullhit~=2,1));   % standard deviation of pulls (1st coord)
    hist(pull(L&pullhit~=2,1),linspace(-limpull-add,limpull+add));
    % draws the histogram
    xlabel(['\mu=',num2str(bpull_mean(1),'%5.3f'),', \sigma=',...
        num2str(bpull_std(1),'%5.3f')]);
    ylims=get(gca,'ylim'); % limits of the plot in y
    hold on
    % plots standard deviation as a red line
    plot([-bpull_std(1) -bpull_std(1)],ylims,'r');
    plot([bpull_std(1) bpull_std(1)],ylims,'r');
    % plots mean as a green line
    plot([bpull_mean(1) bpull_mean(1)],ylims,'g');
end
ylabel('R\Phi');
hold off

% z-pulls at 1st detector
subplot(2,5,2)  % pick upper right corner
box on;
if sum(L)~=0
    % mean and standard deviation
    bpull_mean(2)=mean(pull(L&pullhit~=1,2)); % mean of pulls (2nd coord)
    bpull_std(2)=std(pull(L&pullhit~=1,2));   % standard deviation of pulls (2nd coord)

    hist(pull(L&pullhit~=1,2),linspace(-limpull-add,limpull+add))
    xlabel(['\mu=',num2str(bpull_mean(2),'%5.3f'),', \sigma=',...
        num2str(bpull_std(2),'%5.3f')]);
    ylims=get(gca,'ylim'); % limits of the plot in y
    hold on
    % plots standard deviation as a red line
    plot([-bpull_std(2) -bpull_std(2)],ylims,'r');
    plot([bpull_std(2) bpull_std(2)],ylims,'r');
    % plots mean as a green line
    plot([bpull_mean(2) bpull_mean(2)],ylims,'g');
end
ylabel('z');
hold off

% u-pulls at 1st detector
subplot(2,5,4)  % pick upper left corner
box on;
L=(pullhit>0)&pulltype==0;
if sum(L)~=0
    % mean and standard deviation
    fpull_mean(1)=mean(pull(L&pullhit~=2,1)); % mean of pulls (1st coord)
    fpull_std(1)=std(pull(L&pullhit~=2,1));% std dev of pulls (1st coord)

    hist(pull(L&pullhit~=2,1),linspace(-limpull-add,limpull+add))
    ylabel(['\mu=',num2str(fpull_mean(1),'%5.3f'),', \sigma=',...
        num2str(fpull_std(1),'%5.3f')]);
    ylims=get(gca,'ylim'); % limits of the plot in y
    hold on
    % plots standard deviation as a red line
    plot([-fpull_std(1) -fpull_std(1)],ylims,'r');
    plot([fpull_std(1) fpull_std(1)],ylims,'r');
    % plots mean as a green line
    plot([fpull_mean(1) fpull_mean(1)],ylims,'g');
end
ylabel('u');
hold off

% v-pulls at 1st detector
subplot(2,5,5)  % pick upper right corner
box on;
if sum(L)~=0
    % mean and standard deviation
    fpull_mean(2)=mean(pull(L&pullhit~=1,2)); % mean of pulls (2nd coord)
    fpull_std(2)=std(pull(L&pullhit~=1,2));% std dev of pulls (2nd coord)

    hist(pull(L&pullhit~=1,2),linspace(-limpull-add,limpull+add))
    xlabel(['\mu=',num2str(fpull_mean(2),'%5.3f'),', \sigma=',...
        num2str(fpull_std(2),'%5.3f')]);
    ylims=get(gca,'ylim'); % limits of the plot in y
    hold on
    % plots standard deviation as a red line
    plot([-fpull_std(2) -fpull_std(2)],ylims,'r');
    plot([fpull_std(2) fpull_std(2)],ylims,'r');
    % plots mean as a green line
    plot([fpull_mean(2) fpull_mean(2)],ylims,'g');
end
ylabel('v');
hold off

% MC-pulls at beamtube
L=MCpullhit;
MCpulltitles={'R\Phi','z','\vartheta','\beta','\kappa',};
for k=1:5
    subplot(2,5,k+5)
    % pick five positions in the second row of the figure, one after
    % the other
    hist(real(MC_pull(L,k)),linspace(-limpull-add,limpull+add))
    % draws the histogram
    xlabel(['\mu=',num2str(MCpull_mean(k),'%5.3f'),', \sigma=',...
        num2str(MCpull_std(k),'%5.3f')]);
    ylims=get(gca,'ylim'); % limits of the plot in y
    hold on
    % plots standard deviation as a red line
    plot([-MCpull_std(k) -MCpull_std(k)],ylims,'r');
    plot([MCpull_std(k) MCpull_std(k)],ylims,'r');
    % plots mean as a green line
    plot([MCpull_mean(k) MCpull_mean(k)],ylims,'g');
    hold off
    ylabel(MCpulltitles{k});
end     % for k=1:5 - end

subplot(2,5,3)% pick region in the upper midth to print explaining text
axis off
text(20,0,' ');
h=text(-0.1,0.9,'Top:');
set(h,'fontweight','bold')
text(-0.1,0.8,'Pulls of the ');
text(-0.1,0.75,'innermost hit');
h=text(-0.1,0.35,'Bottom:');
set(h,'fontweight','bold')
text(-0.1,0.25,'Standardized residuals');
text(-0.1,0.2,'of the five track');
text(-0.1,0.15,'parameters at the');
text(-0.1,0.1,'beamtube');