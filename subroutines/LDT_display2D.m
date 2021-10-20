function LDT_display2D(varargin)
%(nameb,namef,VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,SPR,theta1,...
 %   theta2,theta3,theta4,display)
%
% function GEOMETRY3D
% Called by LDT_main
%
% Input:    nameb:  Name of barrel arrangement
%           namef:  Name of forward arrangement
%           VTX:    information of the Vertex Detector
%           SIT:    information of the Silicon Inner Tracker
%           TPC:    information of the Time Projection Chamber
%           SET:    information of the Silicon External Tracker
%           FM1:    information of Forward Module 1
%           FM2:    information of Forward Module 2
%           RM1:    information of Rear Module 1
%           RM2:    information of Rear Module 2
%           SPR:    Start parameter range
%           theta1: theta<theta1: only forward region
%           theta2, theta3: theta2<theta<theta3: only barrel region
%           theta4: theta>theta4: only rear region
%           
% Output:   none
%
% GEOMETRY3D displays a sketch of the the chosen detector arrangement

%global fidlog whandle mhandle disfig hisfig unit Flags
global unit Flags SPR octave
warning off;
close(figure(1));
%close(figure(5));
if octave
    [SPR,N,nameb,namef]=paramconversion;
else
    [SPR,N]=paramconversion;
    nameb=varargin{1};
    namef=varargin{2};
end
%[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=LDT_ReadGeometry(nameb{1},namef{1});
[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,nameb{1},namef{1}]=geomconversion(nameb{1},namef{1});
merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2);
[theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta

drawnames=1;	% set to 0 to avoid drawing of detector names
drawtheta=1;	% set to 0 to avoid drawing of limiting angle's values
if Flags.ScaleDownTPC
    drawwholeTPC=1; % set to 1 to have the program display all TPC layers
else
    drawwholeTPC=0;
end

add1=2;
add2=40;
fsize1=7;
fsize2=8;
dots=100;
dots3D=20;

% find greatest values to display
maxupz1=max(VTX.Length1);maxlowz1=min(VTX.Length2);maxR1=max(VTX.Radius);
maxupz2=max(SIT.Length1);maxlowz2=min(SIT.Length2);maxR2=max(SIT.Radius);
maxupz3=max(TPC.Length1);maxlowz3=min(TPC.Length2);maxR3=max(TPC.Radius);
maxupz4=max(FM1.z);maxlowz4=min(FM1.z);maxR4=max(FM1.outerRadius);
maxupz5=max(FM2.z);maxlowz5=min(FM2.z);maxR5=max(FM2.outerRadius);
maxupz6=max(RM1.z);maxlowz6=min(RM1.z);maxR6=max(RM1.outerRadius);
maxupz7=max(RM2.z);maxlowz7=min(RM2.z);maxR7=max(RM2.outerRadius);
maxupz8=max(SET.Length1);maxlowz8=min(SET.Length2);maxR8=max(SET.Radius);

maxlowz=min([maxlowz1,maxlowz2,maxlowz3,maxlowz4,maxlowz5,maxlowz6,...
    maxlowz7,maxlowz8])*unit-100;
maxupz=max([maxupz1,maxupz2,maxupz3,maxupz4,maxupz5,maxupz6,...
    maxupz7,maxupz7,maxupz8])*unit+100;
maxR=max([maxR1,maxR2,maxR3,maxR4,maxR5,maxR6,maxR7,maxR8])*unit+100;
thetalim1=atan(maxR/maxupz);        % limiting angle
thetalim2=pi+atan(maxR/maxlowz);

h=figure(1);   % sketch is displayed in figure 1
figure(1);
set(h,'Name',['Detector Arrangement: ',nameb{1},'-',namef{1}]);
%set(h,'Units','Pixels','OuterPosition',[30,70,1200,900])    % position,
%width and height
%set(h,'Units','Pixels','OuterPosition',[30,70,900,700])    % position, width and height
%set(h,'DefaultAxesFontSize',9); set(h,'DefaultTextFontSize',9);

% draw axis
x=linspace(maxlowz-0.1,maxupz+0.1,dots*8);
y=linspace(-0.1,maxR+0.1,dots*4);
plot(x,zeros(1,length(x)), 'color','red','linestyle','-.','linewidth',2); % plots the y-axis
%plot(x,0,'r'); % plots the y-axis
if ~octave
    title({'Detector Arrangement:',nameb{1},namef{1}}); % the name of the arrangement
end
hold on; 
%str=['Detector Arrangement: ',nameb{1},'-',namef{1}]; % the name of the arrangement
%text(-100,maxR+100,str);
%text(-500,maxR+500,'Detector Arrangement:');
%text(-500,maxR+300,nameb{1});
%text(-500,maxR+100,namef{1});
                                     % displayed in the title
plot(zeros(1,length(y)),y,'color','red','linestyle','-.','linewidth',2);  % plots the x-axis
%plot(0,y,'r');  % plots the x-axis
xlabel('z-axis [mm]');
ylabel('x- & y-axis [mm]');

% draw Vertex Detector
if VTX.Number~=0
    for k=1:VTX.Number
        x=linspace(VTX.Length2(k),VTX.Length1(k),dots)*unit;
        y=VTX.Radius(k)*ones(1,dots)*unit;
        if (VTX.EffRPhi(k)==0 & (VTX.Effz(k)==-1 | VTX.Effz(k)==0))
            %defines passive scatterer
            plot(x,y,'color','k','linewidth',1);   
            % plots the horizontal lines
        else
            plot(x,y,'color','blue','linewidth',2);
        end
        if drawnames
            text(x(end)+add1,y(end),VTX.Name(k),'fontsize',fsize1);
        end
    end
end

% draw Silicon Inner Tracker
if SIT.Number~=0
    for k=1:SIT.Number
        x=linspace(SIT.Length2(k),SIT.Length1(k),dots)*unit;
        y=SIT.Radius(k)*ones(1,dots)*unit;
        if (SIT.EffRPhi(k)==0 & (SIT.Effz(k)==-1 | SIT.Effz(k)==0))
            %defines passive scatterer
            plot(x,y,'color','k','linewidth',1);
            % plots the horizontal lines
        else
            plot(x,y,'color','cyan','linewidth',2);
        end
        if drawnames
            text(x(end)+add1,y(end),SIT.Name(k),'fontsize',fsize1);
        end
    end
end

% draw Time Projection Chamber
% not drawing every layer, only every tenth
if TPC.Number~=0
    if drawwholeTPC==1
        step=1;
    else
        step=10;
    end
    for k=1:step:TPC.Number
        x=linspace(TPC.Length2(k),TPC.Length1(k),dots)*unit;
        y=TPC.Radius(k)*ones(1,dots)*unit;
        if (TPC.EffRPhi(k)==0 & (TPC.Effz(k)==-1 | TPC.Effz(k)==0))
            %defines passive scatterer
            plot(x,y,'color','k','linewidth',1);
               % plots the horizontal lines
        else
            plot(x,y,'color','green','linewidth',2);
        end
    end
    x=linspace(TPC.Length2(end),TPC.Length1(end),dots)*unit;
    y=TPC.Radius(end)*ones(1,dots)*unit;
    if (TPC.EffRPhi(end)==0 & TPC.Effz(end)==-1)
        %defines passive scatterer
        plot(x,y,'color','k','linewidth',1);
           % plots the outermost layer
    else
        plot(x,y,'color','green','linewidth',2);
    end
end

% draw Silicon External Tracker
if SET.Number~=0
    for k=1:SET.Number
        x=linspace(SET.Length2(k),SET.Length1(k),dots)*unit;
        y=SET.Radius(k)*ones(1,dots)*unit;
        if (SET.EffRPhi(k)==0 & (SET.Effz(k)==-1 | SET.Effz(k)==0))
            %defines passive scatterer
            plot(x,y,'color','k','linewidth',1);   
            % plots the horizontal lines
        else
            plot(x,y,'color','red','linewidth',2);
        end
        if drawnames
            text(x(end)+add1,y(end),SET.Name(k),'fontsize',fsize1);
        end
    end
end


%--------------------------------------------------------------------------

% draw Forward Module 1
if FM1.Number~=0
    for k=1:FM1.Number
        x=FM1.z(k)*ones(1,dots)*unit;
        y=linspace(FM1.innerRadius(k),FM1.outerRadius(k),dots)*unit;
        if (FM1.Effu(k)==0 & (FM1.Effv(k)==-1 | FM1.Effv(k)==0))
            plot(x,y,'color','k','linewidth',1);
        else
            plot(x,y,'color','red','linewidth',2);
        end
        if drawnames
            text(x(end),y(end)+add1,FM1.Name(k),'fontsize',fsize1);
        end
    end
end

% draw Forward Module 2
if FM2.Number~=0
    for k=1:FM2.Number
        x=FM2.z(k)*ones(1,dots)*unit;
        y=linspace(FM2.innerRadius(k),FM2.outerRadius(k),dots)*unit;
        if (FM2.Effu(k)==0 & (FM2.Effv(k)==-1 | FM2.Effv(k)==0))
            plot(x,y,'color','k','linewidth',1);
        else
            plot(x,y,'color','yellow','linewidth',2);
        end
        if drawnames
            text(x(end),y(end)+add1,FM2.Name(k),'fontsize',fsize1);
        end
    end
end

% draw Rear Module 1
if RM1.Number~=0
    for k=1:RM1.Number
        x=RM1.z(k)*ones(1,dots)*unit;
        y=linspace(RM1.innerRadius(k),RM1.outerRadius(k),dots)*unit;
        if (RM1.Effu(k)==0 & (RM1.Effv(k)==-1 | RM1.Effv(k)==0))
            plot(x,y,'color','k','linewidth',1);
        else
            plot(x,y,'color','red','linewidth',2);
        end
        if drawnames
            text(x(end),y(end)+add1,RM1.Name(k),'fontsize',fsize1);
        end
    end
end

% draw Rear Module 2
if RM2.Number~=0
    for k=1:RM2.Number
        x=RM2.z(k)*ones(1,dots)*unit;
        y=linspace(RM2.innerRadius(k),RM2.outerRadius(k),dots)*unit;
        if (RM2.Effu(k)==0 & (RM2.Effv(k)==-1 | RM2.Effv(k)==0))
            plot(x,y,'color','k','linewidth',1);
        else
            plot(x,y,'color','yellow','linewidth',2);
        end
        if drawnames
            text(x(end),y(end)+add1,RM2.Name(k),'fontsize',fsize1);
        end
    end
end

% draw start parameter range of theta
a=2;
for d=1:length(SPR.Thetamin)
    for k=1:2
        switch k
            case 1
                theta=SPR.Thetamin(d)*pi/180;
                z=SPR.z(2);
            case 2
                theta=SPR.Thetamax(d)*pi/180;
                z=SPR.z(1);
        end
        if theta<=thetalim1
            x=linspace(z,maxupz,a*dots);
            y=linspace(0,(maxupz-z)*tan(theta),a*dots);
        elseif (theta>thetalim1 & theta<=thetalim2 & theta~=pi/2)
            x=linspace(z,maxR/tan(theta),a*dots);
            y=linspace(0,maxR,a*dots);
        elseif theta==pi/2
            x=z*ones(1,a*dots);
            y=linspace(0,maxR,a*dots);
        else
            x=linspace(z,maxlowz,a*dots);
            y=linspace(0,(maxlowz-z)*tan(theta),a*dots);
        end
        plot(x,y,'color','red','linestyle','--');
    end
end
% draw limiting angles for forward/barrel/rear region
%add2=40;
for k=1:4
    switch k
        case 1
            theta=theta1;
            %xadd=-600-add;
            angle='\theta_1';
        case 2
            theta=theta2;
            %xadd=-600-add;
            angle='\theta_2';
        case 3
            theta=theta3;
            %xadd=add;
            angle='\theta_3';
        case 4
            theta=theta4;
            %xadd=add;
            angle='\theta_4';
    end
    if theta<=thetalim1
        x=linspace(0,maxupz,a*dots);
        y=linspace(0,(maxupz-z)*tan(theta),a*dots);
    elseif (theta>thetalim1 & theta<=thetalim2 & theta~=pi/2)
        x=linspace(0,maxR/tan(theta),a*dots);
        y=linspace(0,maxR,a*dots);
    elseif theta==pi/2
        x=z*ones(1,a*dots);
        y=linspace(0,maxR,a*dots);
    else
        x=linspace(0,maxlowz,a*dots);
        y=linspace(0,(maxlowz-z)*tan(theta),a*dots);
    end
    plot(x,y,'color','green','linestyle','--');
    if drawtheta
        yadd=add2;
        xadd=add2;
        xtext=x(end)+xadd;
        ytext=y(end)+yadd;
        text(xtext,ytext,[angle,'=',num2str(theta*180/pi,'%3.2f'),'deg'],'fontsize',fsize2);
    end
end
axis equal;
hold off;