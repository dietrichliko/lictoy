function LDT_display3D(varargin)
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
warning off;
global unit Flags octave
warning off
%close(figure(1));
close(figure(5));

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


n=dots3D;
m=dots3D;
phi1=linspace(0,3*pi/2,n);
yb=sin(phi1);
zb=cos(phi1);
Z=repmat(zb',1,m);
Cb=ones(n,m,3);

h=figure(5);   % sketch is displayed in figure 1
figure(5);
set(h,'Name','Detector Arrangement 3D');
%set(h,'Units','Pixels','OuterPosition',[200,100,750,650])    % position, width and height
%set(h,'Units','Pixels','OuterPosition',[50,40,900,700])    % position, width and height
%set(h,'DefaultAxesFontSize',9);
%set(h,'DefaultTextFontSize',9);
hold on
%text(-500,maxR+500,'Detector Arrangement:');
%text(-500,maxR+300,nameb{1});
%text(-500,maxR+100,namef{1});
% draw Vertex Detector
if VTX.Number~=0
    for k=1:VTX.Number
        x=linspace(VTX.Length2(k),VTX.Length1(k),n)*unit;
        y=VTX.Radius(k)*yb*unit;
        z=VTX.Radius(k)*Z*unit;

        if (VTX.EffRPhi(k)==0 & (VTX.Effz(k)==-1 | VTX.Effz(k)==0))
            %defines passive scatterer
            c=Cb*0.3;
        else
            c=Cb;
            c(:,:,1)=0;
            c(:,:,2)=0;
            c(:,:,3)=1;
        end
        %keyboard
        surf(x,y,z,c);
        %surface(x,y,z);

    end
end
axis equal
%view([-30,30])
view(-30,30)
%set(gca,'CameraViewAngleMode','Manual')
if ~octave
    title({'Detector Arrangement:',nameb{1},namef{1}}); % the name of the arrangement
end
if SIT.Number~=0
    for k=1:SIT.Number
        x=linspace(SIT.Length2(k),SIT.Length1(k),n)*unit;
        y=SIT.Radius(k)*yb*unit;
        z=SIT.Radius(k)*Z*unit;

        if (SIT.EffRPhi(k)==0 & (SIT.Effz(k)==-1 | SIT.Effz(k)==0))
            %defines passive scatterer
            c=Cb*0.3;
        else
            c=Cb;
            c(:,:,1)=0;
            c(:,:,2)=0.5;
            c(:,:,3)=1;
        end

        surf(x,y,z,c);
        %surf(x,y,z);

    end
end

if TPC.Number~=0

    if drawwholeTPC==1
        step=1;
    else
        step=10;
    end

    for k=1:step:TPC.Number
        x=linspace(TPC.Length2(k),TPC.Length1(k),n)*unit;
        y=TPC.Radius(k)*yb*unit;
        z=TPC.Radius(k)*Z*unit;

        if (TPC.EffRPhi(k)==0 & (TPC.Effz(k)==-1 | TPC.Effz(k)==0))
            %defines passive scatterer
            c=Cb*0.3;
        else
            c=Cb;
            c(:,:,1)=0;
            c(:,:,2)=1;
            c(:,:,3)=0;
        end

        surf(x,y,z,c);
        %surf(x,y,z);
    end
end

if SET.Number~=0
    for k=1:SET.Number
        x=linspace(SET.Length2(k),SET.Length1(k),n)*unit;
        y=SET.Radius(k)*yb*unit;
        z=SET.Radius(k)*Z*unit;

        if (SET.EffRPhi(k)==0 & (SET.Effz(k)==-1 | SET.Effz(k)==0))
            %defines passive scatterer
            c=Cb*0.3;
        else
            c=Cb;
            c(:,:,1)=1;
            c(:,:,2)=0;
            c(:,:,3)=0;
        end

        surf(x,y,z,c);
        %surf(x,y,z);

    end
end

%--------------------------------------------------------------------------
% draw forward region in 3D

phi2=sort(phi1,'descend');
xp=ones(1,2*n);
y1=sin(phi1); y2=sin(phi2);
z1=cos(phi1); z2=cos(phi2);

if FM1.Number~=0
    for k=1:FM1.Number
        x=FM1.z(k)*xp*unit;
        y=[FM1.innerRadius(k)*y1,FM1.outerRadius(k)*y2]*unit;
        z=[FM1.innerRadius(k)*z1,FM1.outerRadius(k)*z2]*unit;
        if (FM1.Effu(k)==0 & (FM1.Effv(k)==-1 | FM1.Effv(k)==0))
            c=ones(1,3)*0.3;
        else
            c=[1,0,0];
        end

        fill3(x,y,z,c);
        %fill(x,y,c);

    end
end

if FM2.Number~=0
    for k=1:FM2.Number
        x=FM2.z(k)*xp*unit;
        y=[FM2.innerRadius(k)*y1,FM2.outerRadius(k)*y2]*unit;
        z=[FM2.innerRadius(k)*z1,FM2.outerRadius(k)*z2]*unit;
        if (FM2.Effu(k)==0 & (FM2.Effv(k)==-1 | FM2.Effv(k)==0))
            c=ones(1,3)*0.3;
        else
            c=[1,1,0];
        end

        fill3(x,y,z,c);
        %fill(x,y,c);
    end
end

if RM1.Number~=0
    for k=1:RM1.Number
        x=RM1.z(k)*xp*unit;
        y=[RM1.innerRadius(k)*y1,RM1.outerRadius(k)*y2]*unit;
        z=[RM1.innerRadius(k)*z1,RM1.outerRadius(k)*z2]*unit;
        if (RM1.Effu(k)==0 & (RM1.Effv(k)==-1 | RM1.Effv(k)==0))
            c=ones(1,3)*0.3;
        else
            c=[1,0,0];
        end

        fill3(x,y,z,c);
        %fill(x,y,c);
    end
end

if RM2.Number~=0
    for k=1:RM2.Number
        x=RM2.z(k)*xp*unit;
        y=[RM2.innerRadius(k)*y1,RM2.outerRadius(k)*y2]*unit;
        z=[RM2.innerRadius(k)*z1,RM2.outerRadius(k)*z2]*unit;
        if (RM2.Effu(k)==0 & (RM2.Effv(k)==-1 | RM2.Effv(k)==0))
            c=ones(1,3)*0.3;
        else
            c=[1,1,0];
        end

        fill3(x,y,z,c);
        %fill(x,y,c);
    end
end
