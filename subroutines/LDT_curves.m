function LDT_curves(fname)

%global SPR Flags unit
global octave
warning off

den=1.5;
factor=2;
warning off;
if nargin==0
    fname='results.mat';
end
load(fname,'resultsofcurves','nameb','namef','SPR','Flags','unit')
close(figure(6));
close(figure(7));

if exist('resultsofcurves')

    %rms(1)=RPhi rms(2)=z rms(3)=theta rms(4)=phi rms(5)=dpt/pt rms(6)=dpt/pt^2
    %imp(1)=mean(ip2d), imp(2)=std(ip2d), imp(3)=mean(ip3d), imp(4)=std(ip3d)
    %resultsofcurves(1,ICurve,IPoint)=rms(1);
    %resultsofcurves(2,ICurve,IPoint)=rms(2);
    %resultsofcurves(3,ICurve,IPoint)=rms(3);
    %resultsofcurves(4,ICurve,IPoint)=rms(4);
    %resultsofcurves(5,ICurve,IPoint)=rms(5);
    %resultsofcurves(6,ICurve,IPoint)=rms(6);
    %resultsofcurves(7,ICurve,IPoint)=imp(2);
    %resultsofcurves(8,ICurve,IPoint)=imp(3);

    resultsofcurves(1,:,:)=resultsofcurves(1,:,:)*1000*unit; % RPhi in [mum]
    resultsofcurves(2,:,:)=resultsofcurves(2,:,:)*1000*unit; % z in [mum]
    resultsofcurves(3,:,:)=resultsofcurves(3,:,:)*1000;      % theta in [mrad]
    resultsofcurves(4,:,:)=resultsofcurves(4,:,:)*1000;       % phi in [mrad]
    %resultsofcurves(6,:,:)=
    resultsofcurves(7,:,:)=resultsofcurves(7,:,:)*1000*unit; % std(imp2D) in [mum]
    resultsofcurves(8,:,:)=resultsofcurves(8,:,:)*1000*unit; % mean(imp3D) in [mum]

    pos1=[50,70,900,650];
    pos2=[70,50,900,650];

    % Flags.Curves==6: theta and geometry varied, draw in dependence on theta
    % Flags.Curves==5: pt and geometry varied, draw in dependence on theta
    % Flags.Curves==4: both pt and theta varied, draw in dependence on theta
    % Flags.Curves==3: both pt and theta varied, draw in dependence on pt
    % Flags.Curves==2: only pt varied
    % Flags.Curves==1: only theta varied

    if (Flags.Curves==2 | Flags.Curves==3 | Flags.Curves==5)
        % draw in dependence of pt
        pmin=SPR.Ptmin;
        pmax=SPR.Ptmax;
        thetamin=SPR.Thetamin;
        thetamax=SPR.Thetamax;

        if length(pmin)>1
            flagpt=1;
            for k=1:length(pmin)
                p(k)=pmin(k)+(pmax(k)-pmin(k))/2;
            end
        else
            flagpt=0;
        end            

        if Flags.Curves<5
            for k=1:length(thetamin)
                leg{k}=[num2str(thetamin(k)),'-',num2str(thetamax(k)),' deg'];
            end
        else
            if flagpt
                for k=1:length(nameb)
                    barrelfile=nameb{k};
                    forwardfile=namef{k};
                    if octave
                        str=[barrelfile(1:end-4),'-',forwardfile(1:end-4)];
                    else
                        str=[barrelfile(6:end-4),'-',forwardfile(6:end-4)];
                    end
                    leg{k}=[str];
                end
            else
                leg=[];
                for k=1:length(nameb)
                    barrelfile=nameb{k};
                    forwardfile=namef{k};
                    if octave
                      str=[num2str(k),': ',barrelfile(1:end-4),'-',forwardfile(1:end-4)];
                    else
                      str=[num2str(k),': ',barrelfile(6:end-4),'-',forwardfile(6:end-4)];
                    end  
                    while (length(str)>length(leg))&k>1
                        leg(:,end+1)=' ';
                    end
                    while (length(str)<length(leg))&k>1
                        str(:,end+1)=' ';
                    end
                    leg=[leg;str];
                end
                if octave
                  leg=['1-',num2str(length(nameb)),': see parameter file'];
                end
            end
        end

        figure(6) % to display rms values of RPhi, z theta, phi
        figure(6)
        set(gcf,'name',['Reconstructed track parameters (',fname,')'],'position',pos1);
        for k=1:4
            subplot(2,2,k)
            data=squeeze(resultsofcurves(k,:,:))';
            if flagpt
                semilogx(p,data,'-o');
            else
                plot(data,'-o');
            end
            if octave
              if flagpt
                xmin=min(p);xmax=max(p);
              else
                xmin=1; xmax=length(data);
              end
              ymin=min(min(data));
              ymax=max(max(data));
              ymax=ymax+(ymax-ymin)/den;
              axis([xmin,xmax,ymin,ymax]);
            end
            switch k
                case 1
                    %title('rms of R\Phi [\mu m]');
                    ylabel('rms of R\Phi [\mu m]');
                    %semilogx(p,squeeze(resultsofcurves(k,:,:))','-o');
                case 2
                    %title('rms of z [\mu m]');
                    ylabel('rms of z [\mu m]');
                    %semilogx(p,squeeze(resultsofcurves(k,:,:))','-o');
                case 3
                    %title('rms of \theta [mrad]');
                    ylabel('rms of \theta [mrad]');
                    %semilogx(p,squeeze(resultsofcurves(k,:,:))','-o');
                case 4
                    %title('rms of \phi [mrad]');
                    ylabel('rms of \phi [mrad]');
                    %semilogx(p,squeeze(resultsofcurves(k,:,:))','-o');
            end
            
            if octave
              legend(leg);
            else
              legend(leg,'location','southoutside');
            end
            if flagpt
                xlabel('Transverse momentum p_t [GeV/c]');
            end
            grid on
        end

        figure(7) % to display rms values of dpt/pt, dpt/pt2, imp2D, mean imp3D
        figure(7)
        set(gcf,'name',['Momentum resolution and impacts (',fname,')'],'position',pos2);
        for k=1:4
            subplot(2,2,k)
            %loglog(p,squeeze(resultsofcurves(k+4,:,:))','-o');
            data=squeeze(resultsofcurves(k+4,:,:))';
            switch k
                 case 1
                     if flagpt
                         loglog(p,data,'-o');
                     else
                         plot(data,'-o');
                     end
                     if octave
                      if flagpt
                        xmin=min(p);xmax=max(p);
                      else
                        xmin=1; xmax=length(data);
                      end
                      ymin=min(min(data));
                      ymax=max(max(data));
                      ymax=ymax*factor;
                      axis([xmin,xmax,ymin,ymax]);
                    end  
                    %title('rms of \Deltap_t/p_t [1]');
                    ylabel('rms of \Deltap_t/p_t [1]');
                case 2
                    if flagpt
                        loglog(p,data,'-o');
                    else
                        plot(data,'-o');
                    end
                    if octave
                      if flagpt
                        xmin=min(p);xmax=max(p);
                      else
                        xmin=1; xmax=length(data);
                      end
                      ymin=min(min(data));
                      ymax=max(max(data));
                      ymax=ymax*factor;
                      axis([xmin,xmax,ymin,ymax]);
                    end  
                    %title('rms of \Deltap_t/p_t^2 [c/GeV]');
                    ylabel('rms of \Deltap_t/p_t^2 [c/GeV]');
                case 3
                    if flagpt
                        semilogx(p,data,'-o');
                    else
                        plot(data,'-o');
                    end
                    if octave
                      if flagpt
                        xmin=min(p);xmax=max(p);
                      else
                        xmin=1; xmax=length(data);
                      end
                      ymin=min(min(data));
                      ymax=max(max(data));
                      ymax=ymax+(ymax-ymin)/den;
                      axis([xmin,xmax,ymin,ymax]);
                    end  
                    %title('std of projected impact [\mu m]');
                    ylabel('std of projected impact [\mu m]');
                case 4
                    if flagpt
                        semilogx(p,data,'-o');
                    else
                        plot(data,'-o');
                    end
                    if octave
                      if flagpt
                        xmin=min(p);xmax=max(p);
                      else
                        xmin=1; xmax=length(data);
                      end
                      ymin=min(min(data));
                      ymax=max(max(data));
                      ymax=ymax+(ymax-ymin)/den;
                      axis([xmin,xmax,ymin,ymax]);
                    end
                    %title('mean of impact in space [\mu m]');
                    ylabel('mean of impact in space [\mu m]');
            end
            
            if octave
              legend(leg);
            else  
              legend(leg,'location','southoutside');
            end
            if flagpt
                xlabel('Transverse momentum p_t [GeV/c]');
            end
            grid on
        end
    else
        % draw in dependence on theta
        pmin=SPR.Ptmin;
        pmax=SPR.Ptmax;
        thetamin=SPR.Thetamin;
        thetamax=SPR.Thetamax;
        for k=1:length(thetamin)
            theta(k)=thetamin(k)+(thetamax(k)-thetamin(k))/2;
        end

        if Flags.Curves<5
            for k=1:length(pmin)
                leg{k}=[num2str(pmin(k)),'-',num2str(pmax(k)),' GeV/c'];
            end
        else
            for k=1:length(nameb)
                barrelfile=nameb{k};
                forwardfile=namef{k};
                if octave
                    str=[barrelfile(1:end-4),'-',forwardfile(1:end-4)];
                else
                    str=[barrelfile(6:end-4),'-',forwardfile(6:end-4)];
                end
                leg{k}=[str];
            end
        end

        figure(6) % to display rms values of RPhi, z theta, phi
        figure(6)
        set(gcf,'name',['Reconstructed track parameters (',fname,')'],'position',pos1);
        for k=1:4
            subplot(2,2,k)
            %plot(theta,squeeze(resultsofcurves(k,:,:))','-o');
            data=squeeze(resultsofcurves(k,:,:))';
            plot(theta,data,'-o');
            if octave
              xmin=min(theta);xmax=max(theta);
              ymin=min(min(data));
              ymax=max(max(data));
              ymax=ymax+(ymax-ymin)/den;
              axis([xmin,xmax,ymin,ymax]);
            end
            if octave
              legend(leg);
            else
              legend(leg,'location','southoutside');
            end
            grid on
            xlabel('Polar angle \theta [deg]');
            switch k
                case 1
                    %title('rms of R\Phi [\mu m]');
                    ylabel('rms of R\Phi [\mu m]');
                case 2
                    %title('rms of z [\mu m]');
                    ylabel('rms of z [\mu m]');
                case 3
                    %title('rms of \theta [mrad]');
                    ylabel('rms of \theta [mrad]');
                case 4
                    %title('rms of \phi [mrad]');
                    ylabel('rms of \phi [mrad]');
            end
        end

        figure(7) % to display rms values of dpt/pt, dpt/pt2, imp2D, mean imp3D
        figure(7)
        set(gcf,'name',['Momentum resolution and impacts (',fname,')'],'position',pos2);
        for k=1:4
            subplot(2,2,k)
            %plot(theta,squeeze(resultsofcurves(k+4,:,:))','-o');
            data=squeeze(resultsofcurves(k+4,:,:))';
            plot(theta,data,'-o');
            if octave
              xmin=min(theta);xmax=max(theta);
              ymin=min(min(data));
              ymax=max(max(data));
              ymax=ymax+(ymax-ymin)/den;
              axis([xmin,xmax,ymin,ymax]);
            end
            if octave
              legend(leg);
            else
              legend(leg,'location','southoutside');
            end
            %legend(leg,'location','southoutside');
            %legend(leg);
            grid on
            xlabel('Polar angle \theta [deg]');
            switch k
                case 1
                    %title('rms of \Deltap_t/p_t [1]');
                    ylabel('rms of \Deltap_t/p_t [1]');
                case 2
                    %title('rms of \Deltap_t/p_t^2 [c/GeV]');
                    ylabel('rms of \Deltap_t/p_t^2 [c/GeV]');
                case 3
                    %title('std of projected impact [\mu m]');
                    ylabel('std of projected impact [\mu m]');
                case 4
                    %title('mean of impact in space [\mu m]');
                    ylabel('mean of impact in space [\mu m]');
            end
        end
    end
else
    error(['The simulation run "',fname,'" was a single run!']);
end