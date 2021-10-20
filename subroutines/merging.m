function count=merging(VD,IT,TPC,FC,FM1,FM2,RM1,RM2,varargin)

% function MERGING
% Called by LDT_main
% Main program: LDT_main
%
% Input:    VD,IT,TPC,FM1,FM2,RM1,RM2: hold the detector information
%           
% Output:   none, conversation with main program via global variables
%
% MERGING distributes the detector information delivered in the structures
% among sorted arrays for the barrel, forward and rear region. Moreover equal
% radii and z positions are changed and the first non passive layers are
% yielded

%
%global Flags
%global N Radius zlength1 zlength2 alpha bXlen bdetnr1
%global effRPhi effz bdistr sig d
%global fsigu fsigv rsigu rsigv fdu fdv rdu rdv fz rz
%global fXlen finnerradius fouterradius fdistr fdeltau fdeltav feffu feffv
%global rXlen rinnerradius routerradius rdistr rdeltau rdeltav reffu reffv
%global fdetnr1 rdetnr1 unit
%global fidlog whandle mhandle
%global blayername flayername rlayername
%
warning off;
global Flags fidlog whandle mhandle
global N radius zpos delta Xlen detnr1
global eff distr sig d
global delta unit name octave

count1=0;
count2=0;
count3=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Merging of barrel data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% merging of all radii and sorting
Radius=[VD.Radius, IT.Radius, TPC.Radius, FC.Radius];   % 3 barrels
[Radius,I]=sort(Radius);                     % sorts ascending
N.BLayer=length(Radius);                     % number of barrel layers

% merging of all values of the 3 barrels and sorting
bXlen=[VD.Xlen IT.Xlen TPC.Xlen FC.Xlen];                      % Xlens
bXlen=bXlen(I);
sig.RPhi0=[VD.SigmaRPhi IT.SigmaRPhi TPC.SigmaRPhi0 FC.SigmaRPhi];  % measurement errors
sig.RPhi0=sig.RPhi0(I);
sig.RPhi1=[zeros(1,VD.Number) zeros(1,IT.Number) TPC.SigmaRPhi1 zeros(1,FC.Number)];  % measurement errors
sig.RPhi1=sig.RPhi1(I);
sig.CdiffRPhi=[zeros(1,VD.Number) zeros(1,IT.Number) TPC.CdiffRPhi zeros(1,FC.Number)];
sig.CdiffRPhi=sig.CdiffRPhi(I);
sig.z0=[VD.Sigmaz IT.Sigmaz TPC.Sigmaz0 FC.Sigmaz];
sig.z0=sig.z0(I);
sig.z1=[zeros(1,VD.Number) zeros(1,IT.Number) TPC.Sigmaz1 zeros(1,FC.Number)];
sig.z1=sig.z1(I);
sig.Cdiffz=[zeros(1,VD.Number) zeros(1,IT.Number) TPC.Cdiffz zeros(1,FC.Number)];
sig.Cdiffz=sig.Cdiffz(I);
alpha=[VD.Alpha IT.Alpha pi/2*ones(1,TPC.Number) FC.Alpha];     % alphas
alpha=alpha(I);
zlength1=[VD.Length1 IT.Length1 TPC.Length1 FC.Length1];          % upper z limits
zlength1=zlength1(I);
zlength2=[VD.Length2 IT.Length2 TPC.Length2 FC.Length2];          % lower z limits
zlength2=zlength2(I);
bdistr=[VD.Distr IT.Distr 2*ones(1,TPC.Number) FC.Distr];        % error distribution
bdistr=bdistr(I);
d.RPhi=[VD.dRPhi IT.dRPhi NaN*ones(1,TPC.Number) FC.dRPhi];     % strip distances
d.RPhi=d.RPhi(I);
d.z=[VD.dz IT.dz NaN*ones(1,TPC.Number) FC.dz];
d.z=d.z(I);
effRPhi=[VD.EffRPhi IT.EffRPhi TPC.EffRPhi FC.EffRPhi];    % inefficiencies
effRPhi=effRPhi(I);
effz=[VD.Effz IT.Effz TPC.Effz FC.Effz];
effz=effz(I);
blayername=[VD.Name IT.Name TPC.Name FC.Name];
blayername=blayername(I);

f=0;
g=0;
N.bPassive=0;
changepos=1e-6;
% equal radii not possible
changecyl=[];
bdetnr1=1;
for k=1:N.BLayer
    if k~=1
        while sum(Radius(k)==(Radius(1:k-1)))
            Radius(k)=Radius(k)+changepos; % equal radii changed
            g=1;
        end
        if g
            changecyl=[changecyl,blayername(k)];
            g=0;
        end
    end
    if ( effRPhi(k)~=0 | (effz(k)~=-1 & effz(k)~=0) )
        if f~=1
            bdetnr1=k; % finds the first non passive layer
            % bdetnr1 needed to determine the first non passive layer, 
            % where the real pull quantities and residuals can be computed
        end
        f=1;
    else
        N.bPassive=N.bPassive+1;  % counts the number of passive layers
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Merging of forward data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hz=[FM1.z,FM2.z,RM1.z,RM2.z];   % z positions
[hz,I]=sort(hz);                % sorts ascending
hinnerradius=...
    [FM1.innerRadius,FM2.innerRadius,RM1.innerRadius,RM2.innerRadius];
hinnerradius=hinnerradius(I);   % inner radii sorted
houterradius=...
    [FM1.outerRadius,FM2.outerRadius,RM1.outerRadius,RM2.outerRadius];
houterradius=houterradius(I);   % outer radii sorted
heffu=[FM1.Effu,FM2.Effu,RM1.Effu,RM2.Effu];
heffu=heffu(I);                 % efficiency of u sorted
heffv=[FM1.Effv,FM2.Effv,RM1.Effv,RM2.Effv];
heffv=heffv(I);                 % efficiency of v sorted
hfdeltau=[FM1.deltau,FM2.deltau,RM1.deltau,RM2.deltau];
hfdeltau=hfdeltau(I);           % coordinate angle of u sorted
hfdeltav=[FM1.deltav,FM2.deltav,RM1.deltav,RM2.deltav];
hfdeltav=hfdeltav(I);           % coordinate angle of v sorted
hfXlen=[FM1.Xlen,FM2.Xlen,RM1.Xlen,RM2.Xlen];
hfXlen=hfXlen(I);               % radiation length sorted
hfdistr=[FM1.Distr,FM2.Distr,RM1.Distr,RM2.Distr];
hfdistr=hfdistr(I);             % error distribution sorted
hsigu=[FM1.Sigmau,FM2.Sigmau,RM1.Sigmau,RM2.Sigmau];
hsigu=hsigu(I);                 % sigma(u) sorted
hsigv=[FM1.Sigmav,FM2.Sigmav,RM1.Sigmav,RM2.Sigmav];
hsigv=hsigv(I);                 % sigma(v) sorted
hdu=[FM1.du,FM2.du,RM1.du,RM2.du];
hdu=hdu(I);                     % strip distance of u sorted
hdv=[FM1.dv,FM2.dv,RM1.dv,RM2.dv];
hdv=hdv(I);                     % strip distance of v sorted
hlayername=[FM1.Name FM2.Name RM1.Name RM2.Name];
hlayername=hlayername(I);

L=hz>0;             % dividing into two parts, z>0, z<0
fz=[NaN,hz(L)];
rz=[NaN,hz(~L)];
finnerradius=[NaN,hinnerradius(L)];
rinnerradius=[NaN,hinnerradius(~L)];
fouterradius=[NaN,houterradius(L)];
routerradius=[NaN,houterradius(~L)];
feffu=[0,heffu(L)];
reffu=[0,heffu(~L)];
feffv=[-1,heffv(L)];
reffv=[-1,heffv(~L)];
fdeltau=[NaN,hfdeltau(L)];
rdeltau=[NaN,hfdeltau(~L)];
fdeltav=[NaN,hfdeltav(L)];
rdeltav=[NaN,hfdeltav(~L)];
fXlen=[bXlen(1),hfXlen(L)];
rXlen=[bXlen(1),hfXlen(~L)];
fdistr=[NaN,hfdistr(L)];
rdistr=[NaN,hfdistr(~L)];
fsigu=[NaN,hsigu(L)];
rsigu=[NaN,hsigu(~L)];
fsigv=[NaN,hsigv(L)];
rsigv=[NaN,hsigv(~L)];
fdu=[NaN,hdu(L)];
rdu=[NaN,hdu(~L)];
fdv=[NaN,hdv(L)];
rdv=[NaN,hdv(~L)];
flayername=[blayername(1),hlayername(L)];
rlayername=[blayername(1),hlayername(~L)];

N.FLayer=length(fz);    % Number of layers in forward direction
N.RLayer=length(rz);    % rear direction

[rz,I]=sort(rz,'descend');      % rear arrays sorted descending
rinnerradius=rinnerradius(I);
routerradius=routerradius(I);
reffu=reffu(I);
reffv=reffv(I);
rdeltau=rdeltau(I);
rdeltav=rdeltav(I);
rXlen=rXlen(I);
rdistr=rdistr(I);
rsigu=rsigu(I);
rsigv=rsigv(I);
rdu=rdu(I);
rdv=rdv(I);
rlayername=rlayername(I);

if Flags.Smear==0             % uses very small measurement errors
    den=10;              % [micrometer]
    den=den*1e-3/unit;  % transformation to units used internally
    %
    %sig.RPhia=sig.RPhia/den;
    %sig.RPhib=sig.RPhib/den;
    %sig.za=sig.za/den;
    %sig.zb=sig.zb/den;
    %d.RPhi=d.RPhi/den;
    %d.z=d.z/den;
    %fsigu=fsigu/den;
    %rsigu=rsigu/den;
    %fsigv=fsigv/den;
    %rsigv=rsigv/den;
    %fdu=fdu/den;
    %rdu=rdu/den;
    %fdv=fdv/den;
    %rdv=rdv/den;
    %
    sig.RPhia(:)=den;
    sig.RPhib(:)=den;
    sig.za(:)=den;
    sig.zb(:)=den;
    d.RPhi(:)=den;
    d.z(:)=den;
    fsigu(:)=den;
    rsigu(:)=den;
    fsigv(:)=den;
    rsigv(:)=den;
    fdu(:)=den;
    rdu(:)=den;
    fdv(:)=den;
    rdv(:)=den;
    % measurement error 0 impossible due to numerical inaccuracy ->
    % little measurement errors have to be included
    % into the reconstruction and therefore also
    % have to be considered in the simulation, to
    % be statistically consistent
end

f=0;
g=0;
N.fPassive=0;
% equal z positions impossible
changefwd=[];
fdetnr1=1;
for k=1:N.FLayer
    if k~=1
        while sum(fz(k)==fz(1:k-1))
            fz(k)=fz(k)+changepos;   % equal z position changed
            g=1;
        end
        if g
            changefwd=[changefwd,flayername(k)];
            %disp(['Warning: z positions of forward layers are equal! z of layer ',...
            %    char(flayername(k)),' changed by ',num2str(changepos),'.']);
            g=0;
        end
    end
    if ( feffu(k)~=0 | (feffv(k)~=-1 & feffv(k)~=0) )
        if f~=1
            fdetnr1=k; % finds the first non passive layer
            % fdetnr1 needed to determine the first non passive layer, 
            % where the real pull quantities and residuals can be computed
        end
        f=1;
    else
        N.fPassive=N.fPassive+1;  % counts the number of passive layers
    end
end

f=0;
g=0;
N.rPassive=0;
% equal z positions impossible
changerear=[];
rdetnr1=1;
for k=1:N.RLayer
    if k~=1
        while sum(rz(k)==rz(1:k-1))
            rz(k)=rz(k)-changepos;   % equal z changed
            g=1;
        end
        if g
            changerear=[changerear,rlayername(k)];
            %disp(['Warning: z positions of rear layers are equal! z of layer ',...
            %    char(rlayername(k)),' changed by -',num2str(changepos),'.']);
            g=0;
        end
    end
    if ( reffu(k)~=0 | (reffv(k)~=-1 & reffv(k)~=0) )
        if f~=1
            rdetnr1=k; % finds the first non passive layer
            % rdetnr1 needed to determine the first non passive layer,
            % where the real pull quantities and residuals can be computed
        end
        f=1;
    else
        N.rPassive=N.rPassive+1;  % counts the number of passive layers
    end
end

add=0;
if ~isempty([changecyl,changefwd,changerear])
    if octave
        fprintf('\n');
    end
    if ~isempty(changecyl)
        
        s=[];
        if length(changecyl)>1
            for k=1:(length(changecyl)-1)
                s=[s,char(changecyl(k)),', '];
            end
            s=[s,char(changecyl(end))];
        else
            s=char(changecyl);
        end
        s=char(s);
        
        if octave
            if length(varargin)==2
                str=['Warning ',varargin{2},': Radius of following barrel layer(s) changed by ',num2str(changepos),': ',s];
            else
                str=['Warning: Radius of following barrel layer(s) changed by ',num2str(changepos),': ',s];
            end
            count1=fprintf('%s\n',str);
            fflush(stdout);
        else
            str=get(whandle,'String');
            str{2}=['Radius of barrel layer ',s,' changed by ',num2str(changepos),'.'];
            set(whandle,'String',str);
            set(whandle,'foregroundcolor','red')
            if exist('fidlog')
                fprintf(fidlog,'%s\n',str{3});
            end
        end
        
       
        %count1=fprintf('%s\n',['Warning: Radius of layer ',s,' changed by ',num2str(changepos),'.']);
        add=1;
    end
    if ~isempty(changefwd)
        
        s=[];
        if length(changefwd)>1
            for k=1:(length(changefwd)-1)
                s=[s,char(changefwd(k)),', '];
            end
            s=[s,char(changefwd(end))];
        else
            s=char(changefwd);
        end
        s=char(s);
        
        if octave
            if length(varargin)==2
                str=['Warning ',varargin{2},': z-pos. of following forward disk(s) changed by ',num2str(changepos),': ',s];
            else
                str=['Warning: z-pos. of following forward disk(s) changed by ',num2str(changepos),': ',s];
            end
            count2=fprintf('%s\n',str);
            fflush(stdout);
        else
            str=get(whandle,'String');
            str{3}=['z-pos. of fwd disk ',s,' changed by ',num2str(changepos),'.'];
            set(whandle,'String',str);
            set(whandle,'foregroundcolor','red')
            if exist('fidlog')
                fprintf(fidlog,'%s\n',str{3});
            end
        end
       
        %count2=fprintf('%s\n',['Warning: z-pos. of fwd disk ',s,' changed by ',num2str(changepos),'.']);
    end
    if ~isempty(changerear)
        
        s=[];
        if length(changerear)>1
            for k=1:(length(changerear)-1)
                s=[s,char(changerear(k)),', '];
            end
            s=[s,char(changerear(end))];
        else
            s=char(changerear);
        end
        s=char(s);
        
        if octave
            if length(varargin)==2
                str=['Warning ',varargin{2},': z-pos. of following rear disk(s) changed by ',num2str(changepos),': ',s];
            else
                str=['Warning: z-pos. of following rear disk(s) changed by ',num2str(changepos),': ',s];
            end
            count3=fprintf('%s\n',str);
            fflush(stdout);
        else
            str=get(whandle,'String');
            str{4}=['z-pos. of rear disk ',s,' changed by -',num2str(changepos),'.'];
            set(whandle,'String',str);
            set(whandle,'foregroundcolor','red')
            if exist('fidlog')
                fprintf(fidlog,'%s\n',str{4});
            end
        end
    end
end

%-------- Collect data in structures ---------

% strip distances
% d.RPhi and d.z already set in section 'barrel data'
d.ru=rdu;
d.rv=rdv;
d.fu=fdu;
d.fv=fdv;
% Sigmas
% sig.RPhia... already set
sig.fu=fsigu;
sig.fv=fsigv;
sig.ru=rsigu;
sig.rv=rsigv;
% Efficiencies
eff.RPhi=effRPhi;
eff.z=effz;
eff.fu=feffu;
eff.fv=feffv;
eff.ru=reffu;
eff.rv=reffv;
% distributions
distr.b=bdistr;
distr.f=fdistr;
distr.r=rdistr;
% Radii
radius.b=Radius;
radius.fout=fouterradius;
radius.fin=finnerradius;
radius.rout=routerradius;
radius.rin=rinnerradius;
% z positions and dimensions
zpos.bmax=zlength1;
zpos.bmin=zlength2;
zpos.r=rz;
zpos.f=fz;
% Radiation length
Xlen.b=bXlen;
Xlen.f=fXlen;
Xlen.r=rXlen;
% stereo angles
delta.fu=fdeltau;
delta.fv=fdeltav;
delta.ru=rdeltau;
delta.rv=rdeltav;
delta.b=alpha;
% innermost active detectors
detnr1.b=bdetnr1;
detnr1.f=fdetnr1;
detnr1.r=rdetnr1;
% layernames
name.b=blayername;
name.f=flayername;
name.r=rlayername;

count=count1+count2+count3;