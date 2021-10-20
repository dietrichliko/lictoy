function [parameters]=randt(SPR,Flags);

% Function randt generates 3 random parameters (theta,Phi,transverse_momentum)
%
% Also the structure SPR is altered in this function (if the Flags is set) 
% since the traverse momentum is different for each track 

% Theta 
%
%costheta=[0,0.7];
warning off;
usecostheta=1;
Theta=SPR.Theta*pi/180;
if usecostheta
    costheta=sort(cos(Theta));
end
Theta=sort(Theta);

if Theta(1)<pi/2;
    if usecostheta
        theta=acos(costheta(1)+rand(1)*(costheta(2)-costheta(1)));
    else
        theta=Theta(1)+rand(1)*(Theta(2)-Theta(1));
    end
else
    if usecostheta
        theta=acos(costheta(1)-rand(1)*(costheta(2)-costheta(1)));
    else
        theta=Theta(2)-rand(1)*(Theta(2)-Theta(1));
    end
end

% to have 2 symmetric slices of forward/intermediate regions
if ((Theta(1) < pi/2) || (Theta(2) > pi/2)) && Flags.ThetaSymmetry
    if rand(1) > 0.5
        theta=pi-theta;
    end
end
%}
%{
%if SPR.Theta(1)<pi/2;
%    theta=SPR.Theta(1)+rand(1)*(SPR.Theta(2)-SPR.Theta(1));
%else
%    theta=SPR.Theta(2)-rand(1)*(SPR.Theta(2)-SPR.Theta(1));
%end
if Theta(1)<pi/2;
    theta=Theta(1)+rand(1)*(Theta(2)-Theta(1));
else
    theta=Theta(2)-rand(1)*(Theta(2)-Theta(1));
end

% to have 2 symmetric slices of forward/intermediate regions
%if ((SPR.Theta(1) < pi/2) || (SPR.Theta(2) > pi/2)) && Flags.ThetaSymmetry
%    if rand(1) > 0.5
%        theta=pi-theta;
%    end
%end
%}
if ((Theta(1) < pi/2) || (Theta(2) > pi/2)) && Flags.ThetaSymmetry
    if rand(1) > 0.5
        theta=pi-theta;
    end
end

% determines the use of absolute momentum 
if Flags.AbsMomentum == 1
    Pt(1)=SPR.Pt(1)*sin(theta);
    Pt(2)=SPR.Pt(2)*sin(theta); 
else
    Pt=SPR.Pt;
end
Pt=sort(Pt);

transverse_momentum = Pt(1)+rand(1)*(Pt(2)-Pt(1));
phi=sort(SPR.phi*pi/180);
Phi=phi(1)+rand(1)*(phi(2)-phi(1));
%Phi=SPR.phi(1)+rand(1)*(SPR.phi(2)-SPR.phi(1));

parameters=[theta,Phi,transverse_momentum];