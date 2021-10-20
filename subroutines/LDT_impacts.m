function LDT_impacts(fname)
% function IMPACTS
% Called by output
% Main program: LDT_main
%
%	INPUT:	ip2d	array of the projected impact parameters
%		ip2mean	mean value of ip2d
%		ip2dstd standard deviation of ip2d
%		ip3d	array of the 3D impact parameters
%		ip2mean	mean value of ip3d
%		ip2dstd standard deviation of ip3d
%
%   IMPACTS displays the histograms of the impact parameters
%   (projected and 3D)

%(ip2d,ip2d_mean,ip2d_std,ip3d,ip3d_mean,ip3d_std,deltazaxis)
warning off;
close(figure(4));
if nargin==0
    fname='results.mat';
end
load(fname,'AMC_tracknr','R','Avertex','param_start','param_fit','SPR',...
            'Cf_store','Aoffset','pulltype','ndf','chi2','pullhit','resout_coordinate',...
            'Bz','unit','convf','N','GeomVersion','runnumber','Flags');
[ip2d,ip2d_mean,ip2d_std,ip3d,ip3d_mean,ip3d_std,ipz,ipmz,ipsdtz,deltazaxis]=...
    impact_par(AMC_tracknr,R,Avertex,param_fit,Cf_store,N);
drawdz=0;
if drawdz==1
    plotn=4;
else
    plotn=3;
end

global unit
	
    figure(4)
    figure(4)
	set(gcf,'Position',[60,90,675,325]); 
	%set(gcf,'defaulttextfontname','lucidatypewriter');
	%set(gcf,'defaultaxesfontname','lucidatypewriter');
    set(gcf,'Name',['Impact parameters (',fname,')']);

	nbins=50;					% number of bins
	c2mu=1000*unit;				% converting all length units to mu
	ip2d=ip2d*c2mu;
	ip2d_mean=ip2d_mean*c2mu;
	ip2d_std=ip2d_std*c2mu;
	ip3d=ip3d*c2mu;
	ip3d_mean=ip3d_mean*c2mu;
	ip3d_std=ip3d_std*c2mu;
	
%-------------------------- PROJECTED IMPACT --------------------------
%----------------------------------------------------------------------
	lim_max=max(ip2d);
	lim_min=min(ip2d);
	bins=linspace(lim_min,lim_max,nbins);
	
    subplot(1,plotn,1); 
    box on;
	
	hist(ip2d,bins);
    % draws the histogram
    ylabel('Projected impact');
    
	% limits of the plot in y
	ylims=get(gca,'ylim'); 
    hold on;
	
    % plots standard deviation as a red line
    plot([-ip2d_std -ip2d_std],ylims,'r');
    plot([ ip2d_std  ip2d_std],ylims,'r');
    % plots mean as a green line
    plot([ip2d_mean ip2d_mean],ylims,'g');
	xlabel('[\mu m]');
    hold off;
	
% -------------------- 3 DIMENSIONAL IMPACT ---------------------------
%----------------------------------------------------------------------
	% takes only absolute values of ip3d into account
	lim_max=max(ip3d);
	lim_min=min(ip3d);
	bins=linspace(lim_min,lim_max,nbins);
	
	subplot(1,plotn,3); 
    box on;
	
    %n=hist(ip2d,bins);
	hist(ip3d,nbins);
    ylabel('Impact in 3D');
    
	% limits of the plot in y
	ylims=get(gca,'ylim'); 
    hold on;
    % plots mean as a green line
    plot([ip3d_mean ip3d_mean],ylims,'g');
	xlabel('[\mu m]');
    hold off;
	
	

%--------------------------------------------------------------------------
% -- Delta z on the axis (approximation) ----------------------------------
%--------------------------------------------------------------------------
if drawdz==1
    deltazaxis=deltazaxis*1000;
    subplot(1,4,4);
    hold on;
    title('delta(z) on axis (approx.)');
    xlabel('[\mu m]');
    n=hist(deltazaxis,nbins);
    hist(deltazaxis,nbins);
    n=max(n);
    dzmean=mean(deltazaxis);
    dzstd=std(deltazaxis);
    plot(ones(1,100)*dzmean,linspace(0,n),'g');
end
%-----------------------------------------------------------------------
% -- WRITING SOME EXPLAINING TEXT IN THE MIDDLE SECTION OF THE FIGURE --
%-----------------------------------------------------------------------
	x_textpos=0.05;
        
	subplot(1,plotn,2);
	axis off;
    %text(x_textpos+10,0.9,' ');
	h=text(x_textpos,0.9,'Left hand side:');
	set(h,'fontweight','bold');
	text(x_textpos,0.8,['mean:    ',num2str(ip2d_mean)]);
	text(x_textpos,0.7,['std:     ',num2str(ip2d_std)]);
	
	h=text(x_textpos,0.6,'Right hand side:');
	set(h,'fontweight','bold');
	text(x_textpos,0.5,['mean:    ',num2str(ip3d_mean)]);
	text(x_textpos,0.4,['std:     ',num2str(ip3d_std)]);

    if drawdz==1
        h=text(x_textpos,0.3,'Outer right side:');
        set(h,'fontweight','bold');
        text(x_textpos,0.2,['mean:    ',num2str(dzmean)]);
        text(x_textpos,0.1,['std:     ',num2str(dzstd)]);
    else
        text(x_textpos,0.1,' ');
    end
    