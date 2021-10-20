function [theta1,theta2,theta3,theta4]=limitangles;

% function LIMITANGLES
% called by LDT_main
% main program: LDT_main
%
% Input:  none, gets data via global variables
%
% Output: theta1, theta2, theta3, theta4
%
%
% LIMITANGLES calculates the angles to divide the whole theta range of 0-pi into
% the regiones forward (theta<theta1),  forward intermediate
% (theta1<theta<theta2), barrel (theta2<theta<theta3),
% rear intermediate (theta3<theta<theta4) and rear (theta4<theta)

global N zpos radius

if (N.BLayer>1 & N.FLayer>1 & N.RLayer>1)
    % every region equipped
    L=zpos.bmax>0;
    L(1)=0;
    theta1=min(atan(radius.b(L)./zpos.bmax(L)));   % theta<theta1: no b hit
    theta2=max(atan(radius.fout./zpos.f));         % theta>theta2: no forward hit
    theta3=pi-max(atan(radius.rout./abs(zpos.r))); % theta<theta3: no rear hit
    L=zpos.bmin<0;
    L(1)=0;
    theta4=pi-min(atan(radius.b(L)./abs(zpos.bmin(L))));
    % theta>theta4: no b hit
end

if (N.BLayer==1 & N.FLayer>1 & N.RLayer>1)
    % only forward/rear region equipped
    theta2=max(atan(radius.fout./zpos.f));         % theta>theta2: no forward hit
    theta3=pi-max(atan(radius.rout./abs(zpos.r))); % theta<theta3: no rear hit
    theta1=theta2;
    theta4=theta3;
end

if (N.BLayer>1 & N.FLayer>1 & N.RLayer==1)
    % rear region not equipped
    L=zpos.bmax>0;
    L(1)=0;
    theta1=min(atan(Radius(L)./zpos.bmax(L)));   % theta<theta1: no b hit
    theta2=max(atan(radius.fout./zpos.f));         % theta>theta2: no forward hit
    L=zlength2<0;
    L(1)=0;
    theta4=pi-min(atan(Radius(L)./abs(zlength2(L))));
    % theta>theta4: no b hit
    theta3=theta4; % theta<theta3: no rear hit
end

if (N.BLayer>1 & N.FLayer==1 & N.RLayer>1)
    % forward region not equipped
    L=zpos.bmax>0;
    L(1)=0;
    theta1=min(atan(Radius(L)./zpos.bmax(L)));   % theta<theta1: no b hit
    theta2=theta1;         % theta>theta2: no forward hit
    theta3=pi-max(atan(radius.rout./abs(zpos.r))); % theta<theta3: no rear hit
    L=zlength2<0;
    L(1)=0;
    theta4=pi-min(atan(Radius(L)./abs(zlength2(L))));
    % theta>theta4: no b hit
end

if (N.BLayer==1 & N.FLayer>1 & N.RLayer==1)
    % only forward region equipped
    theta2=max(atan(radius.fout./zpos.f));         % theta>theta2: no forward hit
    theta1=theta2;
    theta3=theta2;
    theta4=theta2;
end

if (N.BLayer==1 & N.FLayer==1 & N.RLayer>1)
    % only rear region equipped
    theta3=pi-max(atan(radius.rout./abs(zpos.r))); % theta<theta3: no rear hit
    theta1=theta3;
    theta2=theta3;
    theta4=theta3;
end

if (N.BLayer>1 & N.FLayer==1 & N.RLayer==1)
    % only b region equipped
    L=zpos.bmax>0;
    L(1)=0;
    theta1=min(atan(radius.b(L)./zpos.bmax(L)));   % theta<theta1: no b hit
    theta2=theta1;
    L=zpos.bmin<0;
    L(1)=0;
    theta4=pi-min(atan(radius.b(L)./abs(zpos.bmin(L))));
    % theta>theta4: no b hit
    theta3=theta4;
end