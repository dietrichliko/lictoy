function vertprint(N,Q,VTX,MCparam6,param6,C6,filename)

% function VERTPRINT
% Called by vertexout
% Main program: LDT_main
%
% Input:    N           Number of tracks
%           Q           Array holding the charges of the tracks
%           VTX         Vertex position (x,y,z) [cm]
%           MCparam6    Monte carlo truth of cartesian parameters
%           param6      Fitted cartesian parameters
%           C6          Corresponding 6x6 covarinace matrix
%
% Output:   none
%
% VERTCONV writes the in N specified number of tracks to a text file in the
% DATA HARVESTER's comma seperated variables format. The textfile contains
% the Monte Carlo truth of the vertex, the Monte Carlo truth of the
% parameters, the fitted parameters, the charge and the covariance matrix
% elements. It can be read out with said data harvester and fed into the 
% RAVE/VERTIGO vertex fit toolkit.
% 

%global mhandle filename
warning off;
if isempty(filename)
    filename=['LiC_',num2str(N),'_(',num2str(VTX(1)),',',num2str(VTX(2)),...
        ',',num2str(VTX(3)),').txt'];
end
disp(['Writing file: ',filename]);
%disp(['Writing to file: ',filename]);
fid=fopen(filename,'w');
fprintf(fid,['event: id=',num2str(1),'; run=0; tag="LiC";\n']);
fprintf(fid,'# monte carlo truth\n');
fprintf(fid,['event:simvtx: id=SV1; x=',num2str(VTX(1)),'; y=',...
    num2str(VTX(2)),'; z=',num2str(VTX(3)),'.\n']);
for i=1:N
    %if mod(i,100)==0
    %    disp(['Writing coordinates of track ',num2str(i),'...']);
    %end
    fprintf(fid,['event:simtrack: id=ST',num2str(i),'; pid=1; ',...
        'px=',num2str(MCparam6(i,4)),'; py=',num2str(MCparam6(i,5)),...
        '; pz=',num2str(MCparam6(i,6)),'; q=',num2str(Q(i)),...
        '; vtxid=SV1;','x=',num2str(MCparam6(i,1)),'; y=',...
        num2str(MCparam6(i,2)),'; z=',num2str(MCparam6(i,3)),...
        '; rid=RT',num2str(i),'.\n']);
end
fprintf(fid,'# rec tracks\n');
for i=1:N
    %if mod(i,100)==0
    %    disp(['Writing covariance matrix of track ',num2str(i),'...']);
    %end
    C=squeeze(C6(i,:,:));
    fprintf(fid,['event:track: dxx=',num2str(C(1,1)),'; dxy=',...
        num2str(C(1,2)),'; dxz=',num2str(C(1,3)),'; dxpx=',...
        num2str(C(1,4)),'; dxpy=',num2str(C(1,5)),'; dxpz=',...
        num2str(C(1,6)),'; dyy=',num2str(C(2,2)),'; dyz=',...
        num2str(C(2,3)),'; dypx=',num2str(C(2,4)),'; dypy=',...
        num2str(C(2,5)),'; dypz=',num2str(C(2,6)),'; dzz=',...
        num2str(C(3,3)),'; dzpx=',num2str(C(3,4)),'; dzpy=',...
        num2str(C(3,5)),'; dzpz=',num2str(C(3,6)),'; dpxpx=',...
        num2str(C(4,4)),'; dpxpy=',num2str(C(4,5)),'; dpxpz=',...
        num2str(C(4,6)),'; dpypy=',num2str(C(5,5)),'; dpypz=',...
        num2str(C(5,6)),'; dpzpz=',num2str(C(6,6)),'; id=RT',...
        num2str(i),'; px=',num2str(param6(i,4)),'; py=',...
        num2str(param6(i,5)),'; pz=',num2str(param6(i,6)),'; q=',...
        num2str(Q(i)),'; vtxid=-1','; x=',num2str(param6(i,1)),'; y=',...
        num2str(param6(i,2)),'; z=',num2str(param6(i,3)),'.\n']);
end
disp('finished!');