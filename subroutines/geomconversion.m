function [VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2,nameb,namef]=geomconversion(nameb,namef)

% function LDT_GeomConversion
% Called by LDT_ReadGeometry
% Main function: LDT_main
%
% Input:    Valuesb (contains read out values from the barrel input sheet)
%           Valuesf (contains read out values from the forward input sheet)
% Output:   classes VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2 (contain the data of
%               the Vertex Detector (VTX), Silicon Inner Tracker (SIT), Time
%               Projection Chamber (TPC), Silicon External Tracker (SET), 
%               Forward Module 1 (FM1), Forward Module 2 (FM2), 
%               Rear Module 1 (RM1), Rear Module 2 (RM2)
%
% LDT_GeomConversion fills the values read from the desired input sheet into 
% the classes used by the rest of the program.

global Flags GeomVersion
global unit Bz SPR octave



warning off;
if octave
    if nameb(end)~='m'
      nameb(end)=[];
    end      
    fid=fopen(['geom/',nameb]);
    %keyboard
    d=0;
    while 1
        d=d+1;
        str=fgetl(fid);
        if strcmp(num2str(str),'-1')
            break;
        end
        A=split(str,':');
        S=size(A);
        if S(1)==2;
            Valuesb{d}=A(2,:);
        else
            Valuesb{d}=[];
        end
    end
    fclose(fid);
else
    [Text,Valuesb]=textread(nameb,'%s %s','delimiter',':');
end



if octave
    if namef(end)~='m'
      namef(end)=[];
    end      
    fid=fopen(['geom/',namef]);
    d=0;
    while 1
        d=d+1;
        str=fgetl(fid);
        if strcmp(num2str(str),'-1')
            break;
        end
        A=split(str,':');
        S=size(A);
        if S(1)==2;
            Valuesf{d}=A(2,:);
        else
            Valuesf{d}=[];
        end
    end
    fclose(fid);
else
    [Text,Valuesf]=textread(namef,'%s %s','delimiter',':');
end

        % In 'Text' the text of the input sheet is stored, whereas
        % 'Valuesb' and 'Valuesf' hold the input arguments
        % The function 'textread' uses in that case a ':' as delimiter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Converting values of barrel 1 (VTX) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=3;
GeomVersion.BNo=str2num(Valuesb{k});
GeomVersion.Bfilename=nameb;

k=6;

VTX.Number=str2num(Valuesb{k});
if VTX.Number==0 % use empty arrays if no VTX desired
    VTX.Radius=[];
    VTX.Length1=[];
    VTX.Length2=[];
    VTX.EffRPhi=[];
    VTX.Effz=[];
    VTX.Alpha=[];
    VTX.Xlen=[];
    VTX.Distr=[];
    VTX.SigmaRPhi=[];
    VTX.Sigmaz=[];
    VTX.dRPhi=[];
    VTX.dz=[];
    VTX.Name=[];
    k=24;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    if octave
        VTX.Name=split(Valuesb{k},',');
        S=size(VTX.Name);
        VTX.Name=mat2cell(VTX.Name,ones(1,S(1)),S(2));
        for d=1:length(VTX.Name)
            str=VTX.Name{d};
            str(str==' ')=[];
            VTX.Name{d}=str;
        end
    else
        VTX.Name=strread(char(Valuesb{k}),'%s','delimiter',',');
    end
    %if isempty(VTX.Name)
    if length(VTX.Name)~=VTX.Number
        for d=1:VTX.Number
            name=['VTX-',num2str(d)];
            VTX.Name{d}=name;
        end
        VTX.Name=VTX.Name';
    end
    %else
    %    if length(VTX.Name)~=VTX.Number
    %        next=VTX.Number-length(VTX.Name)+1;
    %        for d=next:VTX.Number
    %            VTX.Name{d}=[];
    %        end
    %    end
    VTX.Name=VTX.Name';
    %end
    k=k+1;
    
    % Radius
    VTX.Radius=str2num(Valuesb{k})/unit;
    switch length(VTX.Radius)
        case VTX.Number
            % do nothing
        case 2
            VTX.Radius=linspace(VTX.Radius(1),VTX.Radius(2),VTX.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['VTX Radius (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;
    
    % Upper z limit
    VTX.Length1=str2num(Valuesb{k})/unit;
    switch length(VTX.Length1)
        case VTX.Number
            % do nothing
        case 1
            VTX.Length1=VTX.Length1*ones(1,VTX.Number);
            % every layer has the same length
        case 2
            VTX.Length1(VTX.Number)=VTX.Length1(2);
            R=VTX.Radius;
            VTX.Length1(2:end-1)=VTX.Length1(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(VTX.Length1(end)-VTX.Length1(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['VTX Upper limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Lower z limit
    VTX.Length2=str2num(Valuesb{k})/unit;
    switch length(VTX.Length2)
        case VTX.Number
            % do nothing
        case 1
            VTX.Length2=VTX.Length2*ones(1,VTX.Number);
            % every layer has the same length
        case 2
            VTX.Length2(VTX.Number)=VTX.Length2(2);
            R=VTX.Radius;
            VTX.Length2(2:end-1)=VTX.Length2(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(VTX.Length2(end)-VTX.Length2(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['VTX Lower limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of RPhi measurements (real numbers between 0 and 1)
    VTX.EffRPhi=str2num(Valuesb{k});
    if any((VTX.EffRPhi<0)|(VTX.EffRPhi>1))
        error(['VTX Efficiency RPhi (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(VTX.EffRPhi)
        case VTX.Number
            % do nothing
        case 1
            VTX.EffRPhi(1:VTX.Number)=VTX.EffRPhi;
            % every layer has the same efficiency
        otherwise
            error(['VTX Efficiency RPhi (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of z measurements (real numbers between 0 and 1)
    VTX.Effz=str2num(Valuesb{k});
    for d=1:length(VTX.Effz)
        eff=VTX.Effz(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['VTX Efficiency z (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(VTX.Effz)
        case VTX.Number
            % do nothing
        case 1
            VTX.Effz(1:VTX.Number)=VTX.Effz;
            % every layer has the same efficiency
        otherwise
            error(['VTX Efficiency z (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    %Lpassive=( VTX.EffRPhi==0 & VTX.Effz==-1 );
    Lactive=( VTX.EffRPhi~=0 | VTX.Effz~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
    
    % Stereo angle alpha (always measuring RPhi, z-measurements can be
    % rotated)
    VTX.Alpha=str2num(Valuesb{k});
    switch length(VTX.Alpha)
        case 1
            VTX.Alpha(Lactive)=VTX.Alpha;
            % every layer has the same angle
        case numact
            VTX.Alpha(Lactive)=VTX.Alpha;
        otherwise
            error(['VTX stereo angle alpha (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    VTX.Alpha(~Lactive)=NaN;
    k=k+1;

    % Radiation length
    VTX.Xlen=str2num(Valuesb{k});
    switch length(VTX.Xlen)
        case VTX.Number
            % do nothing
        case 1
            VTX.Xlen=VTX.Xlen*ones(1,VTX.Number);
            % every layer has the same thickness
        otherwise
            error(['VTX Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=VTX.Xlen<1e-11;
    VTX.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    VTX.Distr=str2num(Valuesb{k});
    switch length(VTX.Distr)
        case 1
            VTX.Distr=VTX.Distr*ones(1,VTX.Number);
            % every layer has the same error distribution
        case numact
            VTX.Distr(Lactive)=VTX.Distr;
        otherwise
            error(['VTX Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    VTX.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(VTX.Distr==1);    % Luni holds the layers with unif. distr. errors
    Lnorm=(VTX.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of RPhi
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        SigmaRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(SigmaRPhi)
            %case VTX.Number
            % do nothing
            case 1
                VTX.SigmaRPhi(Lnorm)=SigmaRPhi;
                % every layer has the same sigma
            case sum(Lnorm)
                VTX.SigmaRPhi(Lnorm)=SigmaRPhi;
            otherwise
                error(['VTX Sigma(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        VTX.dRPhi(Lnorm)=NaN;
        k=k+1;

        % sigma of z
        Sigmaz=str2num(Valuesb{k})*1e-3/unit;
        switch length(Sigmaz)
            %   case VTX.Number
            % do nothing
            case 1
                VTX.Sigmaz(Lnorm)=Sigmaz;
                % every layer has the same sigma
            case sum(Lnorm)
                VTX.Sigmaz(Lnorm)=Sigmaz;
            otherwise
                error(['VTX Sigma(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        VTX.dz(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of RPhi
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        dRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(dRPhi)
            %      case VTX.Number
            % do nothing
            case 1
                VTX.dRPhi(Luni)=dRPhi;
                % every layer has the same strip distance
            case sum(Luni)
                VTX.dRPhi(Luni)=dRPhi;
            otherwise
                error(['VTX d(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        VTX.SigmaRPhi(Luni)=NaN;
        k=k+1;

        % Distance between strips of z
        dz=str2num(Valuesb{k})*1e-3/unit;
        switch length(dz)
            %   case VTX.Number
            % do nothing
            case 1
                VTX.dz(Luni)=dz;
                % every layer has the same strip distance
            case sum(Luni)
                VTX.dz(Luni)=dz;
            otherwise
                error(['VTX d(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        VTX.Sigmaz(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    VTX.SigmaRPhi(~Lactive)=NaN;
    VTX.Sigmaz(~Lactive)=NaN;
    VTX.dRPhi(~Lactive)=NaN;
    VTX.dz(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Converting values of barrel 2 (SIT) %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIT.Number=str2num(Valuesb{k});
if SIT.Number==0 % use empty arrays if no SIT desired
    SIT.Radius=[];
    SIT.Length1=[];
    SIT.Length2=[];
    SIT.EffRPhi=[];
    SIT.Effz=[];
    SIT.Alpha=[];
    SIT.Xlen=[];
    SIT.Distr=[];
    SIT.SigmaRPhi=[];
    SIT.Sigmaz=[];
    SIT.dRPhi=[];
    SIT.dz=[];
    SIT.Name=[];
    k=42;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    
    if octave
        SIT.Name=split(Valuesb{k},',');
        S=size(SIT.Name);
        SIT.Name=mat2cell(SIT.Name,ones(1,S(1)),S(2));
        for d=1:length(SIT.Name)
            str=SIT.Name{d};
            str(str==' ')=[];
            SIT.Name{d}=str;
        end
    else
        SIT.Name=strread(char(Valuesb{k}),'%s','delimiter',',');
    end
    %if isempty(SIT.Name)
    if length(SIT.Name)~=SIT.Number
        for d=1:SIT.Number
            name=['SIT-',num2str(d)];
            SIT.Name{d}=name;
        end
        SIT.Name=SIT.Name';
    end
    %else
    %    if length(SIT.Name)~=SIT.Number
    %        next=SIT.Number-length(SIT.Name)+1;
    %        for d=next:SIT.Number
    %            SIT.Name{d}=[];
    %        end
    %    end
    SIT.Name=SIT.Name';
    %end
    k=k+1;
    
    % Radius
    SIT.Radius=str2num(Valuesb{k})/unit;
    switch length(SIT.Radius)
        case SIT.Number
            % do nothing
        case 2
            SIT.Radius=linspace(SIT.Radius(1),SIT.Radius(2),SIT.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['SIT Radius (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Upper z limit
    SIT.Length1=str2num(Valuesb{k})/unit;
    switch length(SIT.Length1)
        case SIT.Number
            % do nothing
        case 1
            SIT.Length1=SIT.Length1*ones(1,SIT.Number);
        case 2
            SIT.Length1(SIT.Number)=SIT.Length1(2);
            R=SIT.Radius;
            SIT.Length1(2:end-1)=SIT.Length1(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(SIT.Length1(end)-SIT.Length1(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['SIT Upper limit of z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Lower z limit
    SIT.Length2=str2num(Valuesb{k})/unit;
    switch length(SIT.Length2)
        case SIT.Number
            % do nothing
        case 1
            SIT.Length2=SIT.Length2*ones(1,SIT.Number);
        case 2
            SIT.Length2(SIT.Number)=SIT.Length2(2);
            R=SIT.Radius;
            SIT.Length2(2:end-1)=SIT.Length2(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(SIT.Length2(end)-SIT.Length2(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['SIT Lower limit of z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of RPhi measurements (real numbers between 0 and 1)
    SIT.EffRPhi=str2num(Valuesb{k});
    if any((SIT.EffRPhi<0)|(SIT.EffRPhi>1))
        error(['SIT Efficiency RPhi (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(SIT.EffRPhi)
        case SIT.Number
            % do nothing
        case 1
            SIT.EffRPhi(1:SIT.Number)=SIT.EffRPhi;
            % every layer has the same efficiency
        otherwise
            error(['SIT Efficiency RPhi (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of z measurements (real numbers between 0 and 1)
    SIT.Effz=str2num(Valuesb{k});
    for d=1:length(SIT.Effz)
        eff=SIT.Effz(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['SIT Efficiency z (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(SIT.Effz)
        case SIT.Number
            % do nothing
        case 1
            SIT.Effz(1:SIT.Number)=SIT.Effz;
            % every layer has the same efficiency
        otherwise
            error(['SIT Efficiency z (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    Lactive=(SIT.EffRPhi~=0 | SIT.Effz~=-1); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
    
    % Stereo angle alpha
    SIT.Alpha=str2num(Valuesb{k});
    switch length(SIT.Alpha)
        case SIT.Number
            % do nothing
        case 1
            SIT.Alpha(1:SIT.Number)=SIT.Alpha;
            % every layer has the same alpha
        case numact
            SIT.Alpha(Lactive)=SIT.Alpha;
        otherwise
            error(['SIT stereo angle alpha (Line ',num2str(k),'): number of input arguments must be 1 or number of active layers!']);
    end
    SIT.Alpha(~Lactive)=NaN;
    k=k+1;

    % Radiation length
    SIT.Xlen=str2num(Valuesb{k});
    switch length(SIT.Xlen)
        case SIT.Number
            % do nothing
        case 1
            SIT.Xlen=SIT.Xlen*ones(1,SIT.Number);
            % every layer has the same thickness
        otherwise
            error(['SIT Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=SIT.Xlen<1e-11;
    SIT.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    SIT.Distr=str2num(Valuesb{k});
    switch length(SIT.Distr)
        case SIT.Number
            % do nothing
        case 1
            SIT.Distr=SIT.Distr*ones(1,SIT.Number);
            % every layer has the same error distribution
        case numact
            SIT.Distr(Lactive)=SIT.Distr;
        otherwise
            error(['SIT Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    SIT.Distr(~Lactive)=NaN;
    k=k+1;

    % Sigmas and strip distances 
    Luni=(SIT.Distr==1);    % Luni holds the layers with uniformly distributed errors
    Lnorm=(SIT.Distr==0);   % Lnorm holds those with normal distributed errors
    
    % sigma of RPhi
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        SigmaRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(SigmaRPhi)
            case 1
                SIT.SigmaRPhi(Lnorm)=SigmaRPhi;
                % every layer has the same sigma
            case sum(Lnorm)
                SIT.SigmaRPhi(Lnorm)=SigmaRPhi;
            otherwise
                error(['SIT Sigma(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        SIT.dRPhi(Lnorm)=NaN;
        k=k+1;

        % sigma of z
        Sigmaz=str2num(Valuesb{k})*1e-3/unit;
        switch length(Sigmaz)
            case 1
                SIT.Sigmaz(Lnorm)=Sigmaz;
                % every layer has the same sigma
            case sum(Lnorm)
                SIT.Sigmaz(Lnorm)=Sigmaz;
            otherwise
                error(['SIT Sigma(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        SIT.dz(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of RPhi
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        dRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(dRPhi)
           case 1
                SIT.dRPhi(Luni)=dRPhi;
                % every layer has the same strip distance
            case sum(Luni)
                SIT.dRPhi(Luni)=dRPhi;
            otherwise
                error(['SIT d(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        SIT.SigmaRPhi(Luni)=NaN;
        k=k+1;

        % Distance between strips of z
        dz=str2num(Valuesb{k})*1e-3/unit;
        switch length(dz)
            case 1
                SIT.dz(Luni)=dz;
                % every layer has the same strip distance
            case sum(Luni)
                SIT.dz(Luni)=dz;
            otherwise
                error(['SIT d(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        SIT.Sigmaz(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    SIT.SigmaRPhi(~Lactive)=NaN;
    SIT.Sigmaz(~Lactive)=NaN;
    SIT.dRPhi(~Lactive)=NaN;
    SIT.dz(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Converting values of barrel 3 (TPC) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TPC.Number=str2num(Valuesb{k});
if TPC.Number==0    % use empty arrays if no TPC
    TPC.Radius=[];
    TPC.Length1=[];
    TPC.Length2=[];
    TPC.EffRPhi=[];
    TPC.Effz=[];
    TPC.Xlen=[];
    TPC.SigmaRPhi0=[];
    TPC.SigmaRPhi1=[];
    TPC.CdiffRPhi=[];
    TPC.Sigmaz0=[];
    TPC.Sigmaz1=[];
    TPC.Cdiffz=[];
    TPC.Name=[];
    k=58;
else

    if Flags.ScaleDownTPC
        new=round(TPC.Number/5);
        scaledown=TPC.Number/new;
        TPC.Number=new;
    end
    
    for d=1:TPC.Number
        name=['TPC-',num2str(d)];
        TPC.Name{d}=name;
    end

    k=k+1;
    % Radius
    TPC.Radius=str2num(Valuesb{k})/unit;
    switch length(TPC.Radius)
        case TPC.Number
            % do nothing
        case 2
            TPC.Radius=linspace(TPC.Radius(1),TPC.Radius(2),TPC.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['TPC Radius (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Upper z limit
    TPC.Length1=str2num(Valuesb{k})/unit;
    switch length(TPC.Length1)
        case TPC.Number
            % do nothing
        case 1
            TPC.Length1=TPC.Length1*ones(1,TPC.Number);
            % every layer has the same length
        case 2
            TPC.Length1(TPC.Number)=TPC.Length1(2);
            R=TPC.Radius;
            TPC.Length1(2:end-1)=TPC.Length1(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(TPC.Length1(end)-TPC.Length1(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['TPC Upper limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;

    % Lower z limit
    TPC.Length2=str2num(Valuesb{k})/unit;
    switch length(TPC.Length2)
        case TPC.Number
            % do nothing
        case 1
            TPC.Length2=TPC.Length2*ones(1,TPC.Number);
            % every layer has the same length
        case 2
            TPC.Length2(TPC.Number)=TPC.Length2(2);
            R=TPC.Radius;
            TPC.Length2(2:end-1)=TPC.Length2(1)+(R(2:end-1)-R(1))/...
            (R(end)-R(1))*(TPC.Length2(end)-TPC.Length2(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['TPC Lower limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of RPhi measurements (real numbers between 0 and 1)
    TPC.EffRPhi=str2num(Valuesb{k});
    if any((TPC.EffRPhi<0)|(TPC.EffRPhi>1))
        error(['TPC Efficiency RPhi (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(TPC.EffRPhi)
        case TPC.Number
            % do nothing
        case 1
            TPC.EffRPhi(1:TPC.Number)=TPC.EffRPhi;
            % every layer has the same efficiency
        otherwise
            error(['TPC Efficiency RPhi (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of z measurements (real numbers between 0 and 1)
    TPC.Effz=str2num(Valuesb{k});
    if any((TPC.Effz<0)|(TPC.Effz>1))
        error(['TPC Efficiency z (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(TPC.Effz)
        case TPC.Number
            % do nothing
        case 1
            TPC.Effz(1:TPC.Number)=TPC.Effz;
            % every layer has the same efficiency
        otherwise
            error(['TPC Efficiency z (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Radiation length
    TPC.Xlen=str2num(Valuesb{k});
    switch length(TPC.Xlen)
        case TPC.Number
            % do nothing
        case 1
            TPC.Xlen=TPC.Xlen*ones(1,TPC.Number);
            % every layer has the same thickness
        otherwise
            error(['TPC Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=TPC.Xlen<1e-11;
    TPC.Xlen(L)=1e-11;
    k=k+1;
    
    % Sigmas
    % The sigmas are computed according to the following formula:
    % sigma^2 = [ sigma1^2 + sigma2^2 * (abs( z - zmax )) ]
    TPC.SigmaRPhi0=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.SigmaRPhi0)
        case TPC.Number
            % do nothing
        case 1
            TPC.SigmaRPhi0(1:TPC.Number)=TPC.SigmaRPhi0;
            % every layer has the same sigma
        otherwise
            error(['TPC sigma0(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    TPC.SigmaRPhi1=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.SigmaRPhi1)
        case TPC.Number
            % do nothing
        case 1
            TPC.SigmaRPhi1(1:TPC.Number)=TPC.SigmaRPhi1;
            % every layer has the same sigma
        otherwise
            error(['TPC sigma1(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    TPC.CdiffRPhi=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.CdiffRPhi)
        case TPC.Number
            % do nothing
        case 1
            TPC.CdiffRPhi(1:TPC.Number)=TPC.CdiffRPhi;
            % every layer has the same sigma
        otherwise
            error(['TPC Cdiff(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    TPC.Sigmaz0=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.Sigmaz0)
        case TPC.Number
            % do nothing
        case 1
            TPC.Sigmaz0(1:TPC.Number)=TPC.Sigmaz0;
            % every layer has the same sigma
        otherwise
            error(['TPC sigma0(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    TPC.Sigmaz1=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.Sigmaz1)
        case TPC.Number
            % do nothing
        case 1
            TPC.Sigmaz1(1:TPC.Number)=TPC.Sigmaz1;
            % every layer has the same sigma
        otherwise
            error(['TPC sigma1(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    TPC.Cdiffz=str2num(Valuesb{k})*1e-3/unit;
    switch length(TPC.Cdiffz)
        case TPC.Number
            % do nothing
        case 1
            TPC.Cdiffz(1:TPC.Number)=TPC.Cdiffz;
            % every layer has the same sigma
        otherwise
            error(['TPC Cdiff(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    if Flags.ScaleDownTPC
        TPC.Xlen=TPC.Xlen*scaledown;
        TPC.SigmaRPhi0=TPC.SigmaRPhi0/sqrt(scaledown);
        TPC.SigmaRPhi1=TPC.SigmaRPhi1/sqrt(scaledown);
        %TPC.CdiffRPhi=TPC.CdiffRPhi/sqrt(scaledown);
        TPC.Sigmaz0=TPC.Sigmaz0/sqrt(scaledown);
        TPC.Sigmaz1=TPC.Sigmaz1/sqrt(scaledown);
        %TPC.Cdiffz=TPC.Cdiffz/sqrt(scaledown);
    end
    k=k+4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Converting values of Field Cages (SET) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SET.Number=str2num(Valuesb{k});
if SET.Number==0 % use empty arrays if no SET desired
    SET.Radius=[];
    SET.Length1=[];
    SET.Length2=[];
    SET.EffRPhi=[];
    SET.Effz=[];
    SET.Alpha=[];
    SET.Xlen=[];
    SET.Distr=[];
    SET.SigmaRPhi=[];
    SET.Sigmaz=[];
    SET.dRPhi=[];
    SET.dz=[];
    SET.Name=[];
    k=76;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    
    if octave
        SET.Name=split(Valuesb{k},',');
        S=size(SET.Name);
        SET.Name=mat2cell(SET.Name,ones(1,S(1)),S(2));
        for d=1:length(SET.Name)
            str=SET.Name{d};
            str(str==' ')=[];
            SET.Name{d}=str;
        end
    else
        SET.Name=strread(char(Valuesb{k}),'%s','delimiter',',');
    end
    %if isempty(SET.Name)
    if length(SET.Name)~=SET.Number
        for d=1:SET.Number
            name=['SET-',num2str(d)];
            SET.Name{d}=name;
        end
        SET.Name=SET.Name';
    end
    
    %end
    %else
    %    if length(SET.Name)~=SET.Number
    %        next=SET.Number-length(SET.Name)+1;
    %        for d=next:SET.Number
    %            SET.Name{d}=[];
    %        end
    %    end
        SET.Name=SET.Name';
    %end
    k=k+1;
    
    % Radius
    SET.Radius=str2num(Valuesb{k})/unit;
    switch length(SET.Radius)
        case SET.Number
            % do nothing
        case 2
            SET.Radius=linspace(SET.Radius(1),SET.Radius(2),SET.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['SET Radius (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;
    
    % Upper z limit
    SET.Length1=str2num(Valuesb{k})/unit;
    switch length(SET.Length1)
        case SET.Number
            % do nothing
        case 1
            SET.Length1=SET.Length1*ones(1,SET.Number);
            % every layer has the same length
        case 2
            SET.Length1(SET.Number)=SET.Length1(2);
            R=SET.Radius;
            SET.Length1(2:end-1)=SET.Length1(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(SET.Length1(end)-SET.Length1(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['SET Upper limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Lower z limit
    SET.Length2=str2num(Valuesb{k})/unit;
    switch length(SET.Length2)
        case SET.Number
            % do nothing
        case 1
            SET.Length2=SET.Length2*ones(1,SET.Number);
            % every layer has the same length
        case 2
            SET.Length2(SET.Number)=SET.Length2(2);
            R=SET.Radius;
            SET.Length2(2:end-1)=SET.Length2(1)+(R(2:end-1)-R(1))/...
                (R(end)-R(1))*(SET.Length2(end)-SET.Length2(1));
            % the length of each layer is computed from its radius, so that
            % you find it on the surface of a cone
        otherwise
            error(['SET Lower limit in z (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of RPhi measurements (real numbers between 0 and 1)
    SET.EffRPhi=str2num(Valuesb{k});
    if any((SET.EffRPhi<0)|(SET.EffRPhi>1))
        error(['SET Efficiency RPhi (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(SET.EffRPhi)
        case SET.Number
            % do nothing
        case 1
            SET.EffRPhi(1:SET.Number)=SET.EffRPhi;
            % every layer has the same efficiency
        otherwise
            error(['SET Efficiency RPhi (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of z measurements (real numbers between 0 and 1)
    SET.Effz=str2num(Valuesb{k});
    for d=1:length(SET.Effz)
        eff=SET.Effz(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['SET Efficiency z (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(SET.Effz)
        case SET.Number
            % do nothing
        case 1
            SET.Effz(1:SET.Number)=SET.Effz;
            % every layer has the same efficiency
        otherwise
            error(['SET Efficiency z (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    %Lpassive=( SET.EffRPhi==0 & SET.Effz==-1 );
    Lactive=( SET.EffRPhi~=0 | SET.Effz~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
    
    % Stereo angle alpha (always measuring RPhi, z-measurements can be
    % rotated)
    SET.Alpha=str2num(Valuesb{k});
    switch length(SET.Alpha)
        case 1
            SET.Alpha(Lactive)=SET.Alpha;
            % every layer has the same angle
        case numact
            SET.Alpha(Lactive)=SET.Alpha;
        otherwise
            error(['SET stereo angle alpha (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    SET.Alpha(~Lactive)=NaN;
    k=k+1;

    % Radiation length
    SET.Xlen=str2num(Valuesb{k});
    switch length(SET.Xlen)
        case SET.Number
            % do nothing
        case 1
            SET.Xlen=SET.Xlen*ones(1,SET.Number);
            % every layer has the same thickness
        otherwise
            error(['SET Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=SET.Xlen<1e-11;
    SET.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    SET.Distr=str2num(Valuesb{k});
    switch length(SET.Distr)
        case 1
            SET.Distr=SET.Distr*ones(1,SET.Number);
            % every layer has the same error distribution
        case numact
            SET.Distr(Lactive)=SET.Distr;
        otherwise
            error(['SET Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    SET.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(SET.Distr==1);    % Luni holds the layers with unif. distr. errors
    Lnorm=(SET.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of RPhi
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        SigmaRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(SigmaRPhi)
            %case SET.Number
            % do nothing
            case 1
                SET.SigmaRPhi(Lnorm)=SigmaRPhi;
                % every layer has the same sigma
            case sum(Lnorm)
                SET.SigmaRPhi(Lnorm)=SigmaRPhi;
            otherwise
                error(['SET Sigma(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        SET.dRPhi(Lnorm)=NaN;
        k=k+1;

        % sigma of z
        Sigmaz=str2num(Valuesb{k})*1e-3/unit;
        switch length(Sigmaz)
            %   case SET.Number
            % do nothing
            case 1
                SET.Sigmaz(Lnorm)=Sigmaz;
                % every layer has the same sigma
            case sum(Lnorm)
                SET.Sigmaz(Lnorm)=Sigmaz;
            otherwise
                error(['SET Sigma(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        SET.dz(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of RPhi
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        dRPhi=str2num(Valuesb{k})*1e-3/unit;
        switch length(dRPhi)
            %      case SET.Number
            % do nothing
            case 1
                SET.dRPhi(Luni)=dRPhi;
                % every layer has the same strip distance
            case sum(Luni)
                SET.dRPhi(Luni)=dRPhi;
            otherwise
                error(['SET d(RPhi) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        SET.SigmaRPhi(Luni)=NaN;
        k=k+1;

        % Distance between strips of z
        dz=str2num(Valuesb{k})*1e-3/unit;
        switch length(dz)
            %   case SET.Number
            % do nothing
            case 1
                SET.dz(Luni)=dz;
                % every layer has the same strip distance
            case sum(Luni)
                SET.dz(Luni)=dz;
            otherwise
                error(['SET d(z) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        SET.Sigmaz(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    SET.SigmaRPhi(~Lactive)=NaN;
    SET.Sigmaz(~Lactive)=NaN;
    SET.dRPhi(~Lactive)=NaN;
    SET.dz(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Converting magnetic field and beam spot %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solenid magnetic field
Bz=str2num(Valuesb{k});
k=k+1;

% Range in x
SPR.x=sort(str2num(Valuesb{k})/unit);
if SPR.x(1)==SPR.x(2)
    SPR.x=SPR.x+[-eps eps];
end
k=k+1;
% Range in y
SPR.y=sort(str2num(Valuesb{k})/unit);
if SPR.y(1)==SPR.y(2)
    SPR.y=SPR.y+[-eps eps];
end
k=k+1;
% Range in z
SPR.z=sort(str2num(Valuesb{k})/unit);
if SPR.z(1)==SPR.z(2)
    SPR.z=SPR.z+[-eps eps];
end
k=k+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Converting values of Forward Module 1 (FM1) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=3;
GeomVersion.FNo=str2num(Valuesf{k});
GeomVersion.Ffilename=namef;

k=6;
FM1.Number=str2num(Valuesf{k});
if FM1.Number==0 % use empty arrays if no FM1 desired
    FM1.z=[];
    FM1.innerRadius=[];
    FM1.outerRadius=[];
    FM1.Effu=[];
    FM1.Effv=[];
    FM1.deltau=[];
    FM1.deltav=[];
    FM1.Xlen=[];
    FM1.Distr=[];
    FM1.Sigmau=[];
    FM1.Sigmav=[];
    FM1.du=[];
    FM1.dv=[];
    FM1.Name=[];
    k=25;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    
    if octave
        FM1.Name=split(Valuesf{k},',');
        S=size(FM1.Name);
        FM1.Name=mat2cell(FM1.Name,ones(1,S(1)),S(2));
        for d=1:length(FM1.Name)
            str=FM1.Name{d};
            str(str==' ')=[];
            FM1.Name{d}=str;
        end
    else
        FM1.Name=strread(char(Valuesf{k}),'%s','delimiter',',');
    end
    %if isempty(FM1.Name)
    if length(FM1.Name)~=FM1.Number
        for d=1:FM1.Number
            name=['FM1-',num2str(d)];
            FM1.Name{d}=name;
        end
        FM1.Name=FM1.Name';
    end
    %else
    %    if length(FM1.Name)~=FM1.Number
    %        keyboard
    %        next=FM1.Number-length(FM1.Name)+1;
    %        for d=next:FM1.Number
    %            FM1.Name{d}=[];
    %        end
    %    end
        FM1.Name=FM1.Name';
    %end
    k=k+1;
    
    % z positions
    FM1.z=str2num(Valuesf{k})/unit;
    switch length(FM1.z)
        case FM1.Number
            % do nothing
        case 2
            FM1.z=linspace(FM1.z(1),FM1.z(2),FM1.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['FM1 z positions (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Inner radius
    FM1.innerRadius=str2num(Valuesf{k})/unit;
    switch length(FM1.innerRadius)
        case FM1.Number
            % do nothing
        case 1
            FM1.innerRadius=FM1.innerRadius*ones(1,FM1.Number);
            % every layer has the same inner radius
        case 2
            FM1.innerRadius=...
                linspace(FM1.innerRadius(1),FM1.innerRadius(2),FM1.Number);
        otherwise
            error(['FM1 Inner radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Outer radius
    FM1.outerRadius=str2num(Valuesf{k})/unit;
    switch length(FM1.outerRadius)
        case FM1.Number
            % do nothing
        case 1
            FM1.outerRadius=FM1.outerRadius*ones(1,FM1.Number);
            % every layer has the same outer radius
        case 2
            FM1.outerRadius=...
                linspace(FM1.outerRadius(1),FM1.outerRadius(2),FM1.Number);
        otherwise
            error(['FM1 Outer radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of u measurements (real numbers between 0 and 1)
    FM1.Effu=str2num(Valuesf{k});
    if any((FM1.Effu<0)|(FM1.Effu>1))
        error(['FM1 Efficiency u (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(FM1.Effu)
        case FM1.Number
            % do nothing
        case 1
            FM1.Effu(1:FM1.Number)=FM1.Effu;
            % every layer has the same efficiency
        otherwise
            error(['FM1 Efficiency u (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of v measurements (real numbers between 0 and 1, or -1)
    FM1.Effv=str2num(Valuesf{k});
    for d=1:length(FM1.Effv)
        eff=FM1.Effv(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['FM1 Efficiency v (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(FM1.Effv)
        case FM1.Number
            % do nothing
        case 1
            FM1.Effv(1:FM1.Number)=FM1.Effv;
            % every layer has the same efficiency
        otherwise
            error(['FM1 Efficiency v (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    Lactive=( FM1.Effu~=0 | FM1.Effv~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
        
    % Coordinate angle delta1 
    FM1.deltau=str2num(Valuesf{k});
    switch length(FM1.deltau)
        case 1
            FM1.deltau(1:FM1.Number)=FM1.deltau;
            % every layer has the same angle
        case (FM1.Number)
            % do nothing
        case numact
            FM1.deltau(Lactive)=FM1.deltau;
        otherwise
            error(['FM1 Angle 1st coord. (u) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    FM1.deltau(~Lactive)=NaN;
    k=k+1;
    
    % Coordinate angle delta2 
    FM1.deltav=str2num(Valuesf{k});
    switch length(FM1.deltav)
        case 1
            FM1.deltav(1:FM1.Number)=FM1.deltav;
            % every layer has the same angle
        case (FM1.Number)
            % do nothing
        case numact
            FM1.deltav(Lactive)=FM1.deltav;
        otherwise
            error(['FM1 Angle 2nd coord. (v) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    FM1.deltav(~Lactive)=NaN;
    k=k+1;
    
    % Radiation length
    FM1.Xlen=str2num(Valuesf{k});
    switch length(FM1.Xlen)
        case FM1.Number
            % do nothing
        case 1
            FM1.Xlen=FM1.Xlen*ones(1,FM1.Number);
            % every layer has the same thickness
        otherwise
            error(['FM1 Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=FM1.Xlen<1e-11;
    FM1.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    FM1.Distr=str2num(Valuesf{k});
    switch length(FM1.Distr)
        case FM1.Number
            % do nothing
        case 1
            FM1.Distr=FM1.Distr*ones(1,FM1.Number);
            % every layer has the same error distribution
        case numact
            FM1.Distr(Lactive)=FM1.Distr;
        otherwise
            error(['FM1 Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    FM1.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(FM1.Distr==1);    % Luni holds layers with unif. distr. errors
    Lnorm=(FM1.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of u
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        Sigmau=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmau)
            %case FM1.Number
            % do nothing
            case 1
                FM1.Sigmau(Lnorm)=Sigmau;
                % every layer has the same sigma
            case sum(Lnorm)
                FM1.Sigmau(Lnorm)=Sigmau;
            otherwise
                error(['FM1 Sigma(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        FM1.du(Lnorm)=NaN;
        k=k+1;

        % sigma of v
        Sigmav=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmav)
            %   case FM1.Number
            % do nothing
            case 1
                FM1.Sigmav(Lnorm)=Sigmav;
                % every layer has the same sigma
            case sum(Lnorm)
                FM1.Sigmav(Lnorm)=Sigmav;
            otherwise
                error(['FM1 Sigma(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        FM1.dv(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of u
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        du=str2num(Valuesf{k})*1e-3/unit;
        switch length(du)
            %      case FM1.Number
            % do nothing
            case 1
                FM1.du(Luni)=du;
                % every layer has the same strip distance
            case sum(Luni)
                FM1.du(Luni)=du;
            otherwise
                error(['FM1 d(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        FM1.Sigmau(Luni)=NaN;
        k=k+1;

        % Distance between strips of v
        dv=str2num(Valuesf{k})*1e-3/unit;
        switch length(dv)
            %   case FM1.Number
            % do nothing
            case 1
                FM1.dv(Luni)=dv;
                % every layer has the same strip distance
            case sum(Luni)
                FM1.dv(Luni)=dv;
            otherwise
                error(['FM1 d(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        FM1.Sigmav(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    FM1.Sigmau(~Lactive)=NaN;
    FM1.Sigmav(~Lactive)=NaN;
    FM1.du(~Lactive)=NaN;
    FM1.dv(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Converting values of Forward Module 2 (FM2) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FM2.Number=str2num(Valuesf{k});
if FM2.Number==0 % use empty arrays if no FM2 desired
    FM2.z=[];
    FM2.innerRadius=[];
    FM2.outerRadius=[];
    FM2.Effu=[];
    FM2.Effv=[];
    FM2.deltau=[];
    FM2.deltav=[];
    FM2.Xlen=[];
    FM2.Distr=[];
    FM2.Sigmau=[];
    FM2.Sigmav=[];
    FM2.du=[];
    FM2.dv=[];
    FM2.Name=[];
    k=44;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    if octave
        FM2.Name=split(Valuesf{k},',');
        S=size(FM2.Name);
        FM2.Name=mat2cell(FM2.Name,ones(1,S(1)),S(2));
        for d=1:length(FM2.Name)
            str=FM2.Name{d};
            str(str==' ')=[];
            FM2.Name{d}=str;
        end
    else
        FM2.Name=strread(char(Valuesf{k}),'%s','delimiter',',');
    end
    %if isempty(FM2.Name)
    if length(FM2.Name)~=FM2.Number
        for d=1:FM2.Number
            name=['FM2-',num2str(d)];
            FM2.Name{d}=name;
        end
        FM2.Name=FM2.Name';
    end
    %else
    %    if length(FM2.Name)~=FM2.Number
    %        next=FM2.Number-length(FM2.Name)+1;
    %        for d=next:FM2.Number
    %            FM2.Name{d}=[];
    %        end
    %    end
    FM2.Name=FM2.Name';
    %end
    k=k+1;
    
    % z positions
    FM2.z=str2num(Valuesf{k})/unit;
    switch length(FM2.z)
        case FM2.Number
            % do nothing
        case 2
            FM2.z=linspace(FM2.z(1),FM2.z(2),FM2.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['FM2 z positions (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Inner radius
    FM2.innerRadius=str2num(Valuesf{k})/unit;
    switch length(FM2.innerRadius)
        case FM2.Number
            % do nothing
        case 1
            FM2.innerRadius=FM2.innerRadius*ones(1,FM2.Number);
            % every layer has the same inner radius
        case 2
            FM2.innerRadius=...
                linspace(FM2.innerRadius(1),FM2.innerRadius(2),FM2.Number);
        otherwise
            error(['FM2 Inner radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Outer radius
    FM2.outerRadius=str2num(Valuesf{k})/unit;
    switch length(FM2.outerRadius)
        case FM2.Number
            % do nothing
        case 1
            FM2.outerRadius=FM2.outerRadius*ones(1,FM2.Number);
            % every layer has the same outer radius
        case 2
            FM2.outerRadius=...
                linspace(FM2.outerRadius(1),FM2.outerRadius(2),FM2.Number);
        otherwise
            error(['FM2 Outer radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of u measurements (real numbers between 0 and 1)
    FM2.Effu=str2num(Valuesf{k});
    if any((FM2.Effu<0)|(FM2.Effu>1))
        error(['FM2 Efficiency u (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(FM2.Effu)
        case FM2.Number
            % do nothing
        case 1
            FM2.Effu(1:FM2.Number)=FM2.Effu;
            % every layer has the same efficiency
        otherwise
            error(['FM2 Efficiency u (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of v measurements (real numbers between 0 and 1)
    FM2.Effv=str2num(Valuesf{k});
    for d=1:length(FM2.Effv)
        eff=FM2.Effv(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['FM2 Efficiency v (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(FM2.Effv)
        case FM2.Number
            % do nothing
        case 1
            FM2.Effv(1:FM2.Number)=FM2.Effv;
            % every layer has the same efficiency
        otherwise
            error(['FM2 Efficiency v (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    Lactive=( FM2.Effu~=0 | FM2.Effv~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
        
     % Coordinate angle delta1 
    FM2.deltau=str2num(Valuesf{k});
    switch length(FM2.deltau)
        case 1
            FM2.deltau(1:FM2.Number)=FM2.deltau;
            % every layer has the same angle
        case (FM2.Number)
            % do nothing
        case numact
            FM2.deltau(Lactive)=FM2.deltau;
        otherwise
            error(['FM2 Angle 1st coord. (u) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    FM2.deltau(~Lactive)=NaN;
    k=k+1;
    
    % Coordinate angle deltav 
    FM2.deltav=str2num(Valuesf{k});
    switch length(FM2.deltav)
        case 1
            FM2.deltav(1:FM2.Number)=FM2.deltav;
            % every layer has the same angle
        case (FM2.Number)
            % do nothing
        case numact
            FM2.deltav(Lactive)=FM2.deltav;
        otherwise
            error(['FM2 Angle 2nd coord. (v) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    FM2.deltav(~Lactive)=NaN;
    k=k+1;
    
    % Radiation length
    FM2.Xlen=str2num(Valuesf{k});
    switch length(FM2.Xlen)
        case FM2.Number
            % do nothing
        case 1
            FM2.Xlen=FM2.Xlen*ones(1,FM2.Number);
            % every layer has the same thickness
        otherwise
            error(['FM2 Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=FM2.Xlen<1e-11;
    FM2.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    FM2.Distr=str2num(Valuesf{k});
    switch length(FM2.Distr)
        case FM2.Number
            % do nothing
        case 1
            FM2.Distr=FM2.Distr*ones(1,FM2.Number);
            % every layer has the same error distribution
        case numact
            FM2.Distr(Lactive)=FM2.Distr;
        otherwise
            error(['FM2 Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    FM2.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(FM2.Distr==1);    % Luni holds layers with unif. distr. errors
    Lnorm=(FM2.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of u
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        Sigmau=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmau)
            %case FM2.Number
            % do nothing
            case 1
                FM2.Sigmau(Lnorm)=Sigmau;
                % every layer has the same sigma
            case sum(Lnorm)
                FM2.Sigmau(Lnorm)=Sigmau;
            otherwise
                error(['FM2 Sigma(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        FM2.du(Lnorm)=NaN;
        k=k+1;

        % sigma of v
        Sigmav=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmav)
            %   case FM2.Number
            % do nothing
            case 1
                FM2.Sigmav(Lnorm)=Sigmav;
                % every layer has the same sigma
            case sum(Lnorm)
                FM2.Sigmav(Lnorm)=Sigmav;
            otherwise
                error(['FM2 Sigma(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        FM2.dv(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of u
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        du=str2num(Valuesf{k})*1e-3/unit;
        switch length(du)
            %      case FM2.Number
            % do nothing
            case 1
                FM2.du(Luni)=du;
                % every layer has the same strip distance
            case sum(Luni)
                FM2.du(Luni)=du;
            otherwise
                error(['FM2 d(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        FM2.Sigmau(Luni)=NaN;
        k=k+1;

        % Distance between strips of v
        dv=str2num(Valuesf{k})*1e-3/unit;
        switch length(dv)
            %   case FM2.Number
            % do nothing
            case 1
                FM2.dv(Luni)=dv;
                % every layer has the same strip distance
            case sum(Luni)
                FM2.dv(Luni)=dv;
            otherwise
                error(['FM2 d(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        FM2.Sigmav(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    FM2.Sigmau(~Lactive)=NaN;
    FM2.Sigmav(~Lactive)=NaN;
    FM2.du(~Lactive)=NaN;
    FM2.dv(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Converting values of Rear Module 1 (RM1) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RM1.Number=str2num(Valuesf{k});
if RM1.Number==0 % use empty arrays if no RM1 desired
    RM1.z=[];
    RM1.innerRadius=[];
    RM1.outerRadius=[];
    RM1.Effu=[];
    RM1.Effv=[];
    RM1.deltau=[];
    RM1.deltav=[];
    RM1.Xlen=[];
    RM1.Distr=[];
    RM1.Sigmau=[];
    RM1.Sigmav=[];
    RM1.du=[];
    RM1.dv=[];
    RM1.Name=[];
    k=63;
elseif RM1.Number==-1
    RM1.Number=FM1.Number;
    RM1.z=-FM1.z;
    RM1.innerRadius=FM1.innerRadius;
    RM1.outerRadius=FM1.outerRadius;
    RM1.Effu=FM1.Effu;
    RM1.Effv=FM1.Effv;
    RM1.deltau=FM1.deltau;
    RM1.deltav=FM1.deltav;
    RM1.Xlen=FM1.Xlen;
    RM1.Distr=FM1.Distr;
    RM1.Sigmau=FM1.Sigmau;
    RM1.Sigmav=FM1.Sigmav;
    RM1.du=FM1.du;
    RM1.dv=FM1.dv;
    RM1.Name=FM1.Name;
    for d=1:RM1.Number
       if RM1.Name{d}(1)=='F'
          RM1.Name{d}(1)='R';
       end
    end
    k=63;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    %
    if octave
        RM1.Name=split(Valuesf{k},',');
        S=size(RM1.Name);
        RM1.Name=mat2cell(RM1.Name,ones(1,S(1)),S(2));
        for d=1:length(RM1.Name)
            str=RM1.Name{d};
            str(str==' ')=[];
            RM1.Name{d}=str;
        end
    else
        RM1.Name=strread(char(Valuesf{k}),'%s','delimiter',',');
    end
    %if isempty(RM1.Name)
    if length(RM1.Name)~=RM1.Number
        for d=1:RM1.Number
            name=['RM1-',num2str(d)];
            RM1.Name{d}=name;
        end
        RM1.Name=RM1.Name';
    end
    %else
    %    if length(RM1.Name)~=RM1.Number
    %        next=RM1.Number-length(RM1.Name)+1;
    %        for d=next:RM1.Number
    %            RM1.Name{d}=[];
    %        end
    %    end
    RM1.Name=RM1.Name';
    %end
    k=k+1;
    
    % z positions
    RM1.z=str2num(Valuesf{k})/unit;
    switch length(RM1.z)
        case RM1.Number
            % do nothing
        case 2
            RM1.z=linspace(RM1.z(1),RM1.z(2),RM1.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['RM1 z positions (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Inner radius
    RM1.innerRadius=str2num(Valuesf{k})/unit;
    switch length(RM1.innerRadius)
        case RM1.Number
            % do nothing
        case 1
            RM1.innerRadius=RM1.innerRadius*ones(1,RM1.Number);
            % every layer has the same inner radius
        case 2
            RM1.innerRadius=...
                linspace(RM1.innerRadius(1),RM1.innerRadius(2),RM1.Number);
        otherwise
            error(['RM1 Inner radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Outer radius
    RM1.outerRadius=str2num(Valuesf{k})/unit;
    switch length(RM1.outerRadius)
        case RM1.Number
            % do nothing
        case 1
            RM1.outerRadius=RM1.outerRadius*ones(1,RM1.Number);
            % every layer has the same outer radius
        case 2
            RM1.outerRadius=...
                linspace(RM1.outerRadius(1),RM1.outerRadius(2),RM1.Number);
        otherwise
            error(['RM1 Outer radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of u measurements (real numbers between 0 and 1)
    RM1.Effu=str2num(Valuesf{k});
    if any((RM1.Effu<0)|(RM1.Effu>1))
        error(['RM1 Efficiency u (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(RM1.Effu)
        case RM1.Number
            % do nothing
        case 1
            RM1.Effu(1:RM1.Number)=RM1.Effu;
            % every layer has the same efficiency
        otherwise
            error(['RM1 Efficiency u (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of v measurements (real numbers between 0 and 1)
    RM1.Effv=str2num(Valuesf{k});
    for d=1:length(RM1.Effv)
        eff=RM1.Effv(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['RM1 Efficiency v (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(RM1.Effv)
        case RM1.Number
            % do nothing
        case 1
            RM1.Effv(1:RM1.Number)=RM1.Effv;
            % every layer has the same efficiency
        otherwise
            error(['RM1 Efficiency v (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    Lactive=( RM1.Effu~=0 | RM1.Effv~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
        
     % Coordinate angle deltau 
    RM1.deltau=str2num(Valuesf{k});
    switch length(RM1.deltau)
        case 1
            RM1.deltau(1:RM1.Number)=RM1.deltau;
            % every layer has the same angle
        case (RM1.Number)
            % do nothing
        case numact
            RM1.deltau(Lactive)=RM1.deltau;
        otherwise
            error(['RM1 Angle 1st coord. (u) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    RM1.deltau(~Lactive)=NaN;
    k=k+1;
    
    % Coordinate angle deltav 
    RM1.deltav=str2num(Valuesf{k});
    switch length(RM1.deltav)
        case 1
            RM1.deltav(1:RM1.Number)=RM1.deltav;
            % every layer has the same angle
        case (RM1.Number)
            % do nothing
        case numact
            RM1.deltav(Lactive)=RM1.deltav;
        otherwise
            error(['RM1 Angle 2nd coord. (v) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    k=k+1;
    RM1.deltav(~Lactive)=NaN;
    
    % Radiation length
    RM1.Xlen=str2num(Valuesf{k});
    switch length(RM1.Xlen)
        case RM1.Number
            % do nothing
        case 1
            RM1.Xlen=RM1.Xlen*ones(1,RM1.Number);
            % every layer has the same thickness
        otherwise
            error(['RM1 Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=RM1.Xlen<1e-11;
    RM1.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    RM1.Distr=str2num(Valuesf{k});
    switch length(RM1.Distr)
        case RM1.Number
            % do nothing
        case 1
            RM1.Distr=RM1.Distr*ones(1,RM1.Number);
            % every layer has the same error distribution
        case numact
            RM1.Distr(Lactive)=RM1.Distr;
        otherwise
            error(['RM1 Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    RM1.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(RM1.Distr==1);    % Luni holds layers with unif. distr. errors
    Lnorm=(RM1.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of u
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        Sigmau=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmau)
            %case RM1.Number
            % do nothing
            case 1
                RM1.Sigmau(Lnorm)=Sigmau;
                % every layer has the same sigma
            case sum(Lnorm)
                RM1.Sigmau(Lnorm)=Sigmau;
            otherwise
                error(['RM1 Sigma(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        RM1.du(Lnorm)=NaN;
        k=k+1;

        % sigma of v
        Sigmav=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmav)
            %   case RM1.Number
            % do nothing
            case 1
                RM1.Sigmav(Lnorm)=Sigmav;
                % every layer has the same sigma
            case sum(Lnorm)
                RM1.Sigmav(Lnorm)=Sigmav;
            otherwise
                error(['RM1 Sigma(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        RM1.dv(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of u
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        du=str2num(Valuesf{k})*1e-3/unit;
        switch length(du)
            %      case RM1.Number
            % do nothing
            case 1
                RM1.du(Luni)=du;
                % every layer has the same strip distance
            case sum(Luni)
                RM1.du(Luni)=du;
            otherwise
                error(['RM1 d(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        RM1.Sigmau(Luni)=NaN;
        k=k+1;

        % Distance between strips of v
        dv=str2num(Valuesf{k})*1e-3/unit;
        switch length(dv)
            %   case RM1.Number
            % do nothing
            case 1
                RM1.dv(Luni)=dv;
                % every layer has the same strip distance
            case sum(Luni)
                RM1.dv(Luni)=dv;
            otherwise
                error(['RM1 d(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        RM1.Sigmav(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    RM1.Sigmau(~Lactive)=NaN;
    RM1.Sigmav(~Lactive)=NaN;
    RM1.du(~Lactive)=NaN;
    RM1.dv(~Lactive)=NaN;
    k=k+3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Converting values of Rear Module 2 (RM2) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RM2.Number=str2num(Valuesf{k});
if RM2.Number==0 % use empty arrays if no RM2 desired
    RM2.z=[];
    RM2.innerRadius=[];
    RM2.outerRadius=[];
    RM2.Effu=[];
    RM2.Effv=[];
    RM2.deltau=[];
    RM2.deltav=[];
    RM2.Xlen=[];
    RM2.Distr=[];
    RM2.Sigmau=[];
    RM2.Sigmav=[];
    RM2.du=[];
    RM2.dv=[];
    RM2.Name=[];
    k=82;
elseif RM2.Number==-1
    RM2.Number=FM2.Number;
    RM2.z=-FM2.z;
    RM2.innerRadius=FM2.innerRadius;
    RM2.outerRadius=FM2.outerRadius;
    RM2.Effu=FM2.Effu;
    RM2.Effv=FM2.Effv;
    RM2.deltau=FM2.deltau;
    RM2.deltav=FM2.deltav;
    RM2.Xlen=FM2.Xlen;
    RM2.Distr=FM2.Distr;
    RM2.Sigmau=FM2.Sigmau;
    RM2.Sigmav=FM2.Sigmav;
    RM2.du=FM2.du;
    RM2.dv=FM2.dv;
    RM2.Name=FM2.Name;
    for d=1:RM2.Number
           if RM2.Name{d}(1)=='F'
              RM2.Name{d}(1)='R';
           end
    end
    k=82;
else
    k=k+1;
    
    % Description (optional)
    k=k+1;
    
    % Names of the layers (opt.)
    %
    if octave
        RM2.Name=split(Valuesf{k},',');
        S=size(RM2.Name);
        RM2.Name=mat2cell(RM2.Name,ones(1,S(1)),S(2));
        for d=1:length(RM2.Name)
            str=RM2.Name{d};
            str(str==' ')=[];
            RM2.Name{d}=str;
        end
    else
        RM2.Name=strread(char(Valuesf{k}),'%s','delimiter',',');
    end
    %if isempty(RM2.Name)
    if length(RM2.Name)~=RM2.Number
        for d=1:RM2.Number
            name=['RM2-',num2str(d)];
            RM2.Name{d}=name;
        end
        RM2.Name=RM2.Name';
    end
    %else
    %    if length(RM2.Name)~=RM2.Number
    %        next=RM2.Number-length(RM2.Name)+1;
    %        for d=next:RM2.Number
    %            RM2.Name{d}=[];
    %        end
    %    end
    RM2.Name=RM2.Name';
    %end
    k=k+1;
    
    % z positions
    RM2.z=str2num(Valuesf{k})/unit;
    switch length(RM2.z)
        case RM2.Number
            % do nothing
        case 2
            RM2.z=linspace(RM2.z(1),RM2.z(2),RM2.Number);
            % The interval between the lower and the upper input argument
            % is devided evenly among the desired number of layers
        otherwise
            error(['RM2 z positions (Line ',num2str(k),'): number of input arguments has to be 2 or number of layers!']);
    end
    k=k+1;

    % Inner radius
    RM2.innerRadius=str2num(Valuesf{k})/unit;
    switch length(RM2.innerRadius)
        case RM2.Number
            % do nothing
        case 1
            RM2.innerRadius=RM2.innerRadius*ones(1,RM2.Number);
            % every layer has the same inner radius
        case 2
            RM2.innerRadius=...
                linspace(RM2.innerRadius(1),RM2.innerRadius(2),RM2.Number);
        otherwise
            error(['RM2 Inner radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Outer radius
    RM2.outerRadius=str2num(Valuesf{k})/unit;
    switch length(RM2.outerRadius)
        case RM2.Number
            % do nothing
        case 1
            RM2.outerRadius=RM2.outerRadius*ones(1,RM2.Number);
            % every layer has the same outer radius
        case 2
            RM2.outerRadius=...
                linspace(RM2.outerRadius(1),RM2.outerRadius(2),RM2.Number);
        otherwise
            error(['RM2 Outer radius (Line ',num2str(k),'): number of input arguments has to be 1 or 2 or number of layers!']);
    end
    k=k+1;
    
    % Efficiency of u measurements (real numbers between 0 and 1)
    RM2.Effu=str2num(Valuesf{k});
    if any((RM2.Effu<0)|(RM2.Effu>1))
        error(['RM2 Efficiency u (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1!']);
    end
    switch length(RM2.Effu)
        case RM2.Number
            % do nothing
        case 1
            RM2.Effu(1:RM2.Number)=RM2.Effu;
            % every layer has the same efficiency
        otherwise
            error(['RM2 Efficiency u (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
        
    % Efficiency of v measurements (real numbers between 0 and 1)
    RM2.Effv=str2num(Valuesf{k});
    for d=1:length(RM2.Effv)
        eff=RM2.Effv(d);
        if ( (eff<0 & eff~=-1) | (eff>1) )
            error(['RM2 Efficiency v (Line ',num2str(k),'): input arguments have to be real numbers between 0 and 1 or -1!']);
        end
    end
    switch length(RM2.Effv)
        case RM2.Number
            % do nothing
        case 1
            RM2.Effv(1:RM2.Number)=RM2.Effv;
            % every layer has the same efficiency
        otherwise
            error(['RM2 Efficiency v (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    k=k+1;
    
    Lactive=( RM2.Effu~=0 | RM2.Effv~=-1 ); % 0 if passive, 1 if active
    numact=sum(Lactive);                      % number of active layers
        
     % Coordinate angle deltau 
    RM2.deltau=str2num(Valuesf{k});
    switch length(RM2.deltau)
        case 1
            RM2.deltau(1:RM2.Number)=RM2.deltau;
            % every layer has the same angle
        case (RM2.Number)
            % do nothing
        case numact
            RM2.deltau(Lactive)=RM2.deltau;
        otherwise
            error(['RM2 Angle 1st coord. (u) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    RM2.deltau(~Lactive)=NaN;
    k=k+1;
    
    % Coordinate angle deltav 
    RM2.deltav=str2num(Valuesf{k});
    switch length(RM2.deltav)
        case 1
            RM2.deltav(1:RM2.Number)=RM2.deltav;
            % every layer has the same angle
        case (RM2.Number)
            % do nothing
        case numact
            RM2.deltav(Lactive)=RM2.deltav;
        otherwise
            error(['RM2 Angle 2nd coord. (v) (Line ',num2str(k),'): number of input arguments must be 1 or equal the number of active layers!']);
    end
    RM2.deltav(~Lactive)=NaN;
    k=k+1;
    
    % Radiation length
    RM2.Xlen=str2num(Valuesf{k});
    switch length(RM2.Xlen)
        case RM2.Number
            % do nothing
        case 1
            RM2.Xlen=RM2.Xlen*ones(1,RM2.Number);
            % every layer has the same thickness
        otherwise
            error(['RM2 Thickness (Line ',num2str(k),'): number of input arguments has to be 1 or number of layers!']);
    end
    % due to numerical inaccuracies no Xlen<1e-11 is possible 
    L=RM2.Xlen<1e-11;
    RM2.Xlen(L)=1e-11;
    k=k+1;

    % Error distribution (0: normal, 1: uniform)
    RM2.Distr=str2num(Valuesf{k});
    switch length(RM2.Distr)
        case RM2.Number
            % do nothing
        case 1
            RM2.Distr=RM2.Distr*ones(1,RM2.Number);
            % every layer has the same error distribution
        case numact
            RM2.Distr(Lactive)=RM2.Distr;
        otherwise
            error(['RM2 Error distribution (Line ',num2str(k),'): number of input arguments has to be 1 or number of active layers!']);
    end
    RM2.Distr(~Lactive)=NaN;
    k=k+1;
    
    % Sigmas and strip distances 
    Luni=(RM2.Distr==1);    % Luni holds layers with unif. distr. errors
    Lnorm=(RM2.Distr==0);   % Lnorm holds those with normal distr. errors
    
    % sigma of u
    if sum(Luni)~=numact    % skip this part if only uniform errors desired 
        Sigmau=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmau)
            %case RM2.Number
            % do nothing
            case 1
                RM2.Sigmau(Lnorm)=Sigmau;
                % every layer has the same sigma
            case sum(Lnorm)
                RM2.Sigmau(Lnorm)=Sigmau;
            otherwise
                error(['RM2 Sigma(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        RM2.du(Lnorm)=NaN;
        k=k+1;

        % sigma of v
        Sigmav=str2num(Valuesf{k})*1e-3/unit;
        switch length(Sigmav)
            %   case RM2.Number
            % do nothing
            case 1
                RM2.Sigmav(Lnorm)=Sigmav;
                % every layer has the same sigma
            case sum(Lnorm)
                RM2.Sigmav(Lnorm)=Sigmav;
            otherwise
                error(['RM2 Sigma(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active normal layers!']);
        end
        RM2.dv(Lnorm)=NaN;
        k=k+1;
    else
        k=k+2;
    end

    % Distance between strips of u
    if sum(Lnorm)~=numact  % skips this part if only normal errors desired
        du=str2num(Valuesf{k})*1e-3/unit;
        switch length(du)
            %      case RM2.Number
            % do nothing
            case 1
                RM2.du(Luni)=du;
                % every layer has the same strip distance
            case sum(Luni)
                RM2.du(Luni)=du;
            otherwise
                error(['RM2 d(u) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        RM2.Sigmau(Luni)=NaN;
        k=k+1;

        % Distance between strips of v
        dv=str2num(Valuesf{k})*1e-3/unit;
        switch length(dv)
            %   case RM2.Number
            % do nothing
            case 1
                RM2.dv(Luni)=dv;
                % every layer has the same strip distance
            case sum(Luni)
                RM2.dv(Luni)=dv;
            otherwise
                error(['RM2 d(v) (Line ',num2str(k),'): number of input arguments has to be 1 or number of active uniform layers!']);
        end
        RM2.Sigmav(Luni)=NaN;
        k=k+1;
    else
        k=k+2;
    end
    
    RM2.Sigmau(~Lactive)=NaN;
    RM2.Sigmav(~Lactive)=NaN;
    RM2.du(~Lactive)=NaN;
    RM2.dv(~Lactive)=NaN;
    k=k+3;
end