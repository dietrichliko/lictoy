function [paramf,ierr]=pcyl(R,params,sinbmx)

%		funciton pcyl propagates the state vector from an arbirary vertex 
%		within the beam to the surface of the beam tube


global Bz unit hwait convf

%params=[0,2,0,pi/2,pi/2,0.0001];

ierr=0;
paramf=zeros(1,5);
x_s=params(1);
y_s=params(2);

% setting x_s and y_s to zero for numerical reasons
if abs(x_s)<2*eps
    x_s=0;
end
if abs(y_s)<2*eps
    y_s=0;
end

% error if vertex lies outside of the beamtube
if sqrt(x_s^2+y_s^2)>=R
    ierr=1;
else
    z_s=params(3);
    theta_s=params(4);
    phi_s=params(5);
    kappa_s=params(6);

    h=sign(kappa_s*Bz);            % helicity
    Rh=1/abs(kappa_s);             % Radius of Helix
    Phi_s=atan2(y_s,x_s);          % Polar angle of start point
	
	%here is something wrong!!!
    x_h=x_s-h*Rh*sin(phi_s);       % x-coordinate of helix axis
    y_h=y_s+h*Rh*cos(phi_s);       % y-coordinate of helix axis

    % Coefficients of potence line Ax+By=C
    % The intersection points of the helix and the cylinder are on this line
    A=2*x_h;
    B=2*y_h;
    C=x_h^2+y_h^2-Rh^2+R^2;

    % calculate coordinates of the two intersection points
    if A==0
        y(1:2)=C/B;
        x(1)=sqrt(R^2-y(1).^2);  % y(1)=>y.
        x(2)=-x(1);
    elseif B==0
        x(1:2)=C/A;
        y(1)=sqrt(R^2-x(1).^2);
        y(2)=-y(1);
    elseif abs(A)<=abs(B)
        x(1)=1/(A^2+B^2)*(A*C+B*sqrt((A^2+B^2)*R^2-C^2));
        x(2)=1/(A^2+B^2)*(A*C-B*sqrt((A^2+B^2)*R^2-C^2));
        y=-A/B*x+C/B;
    else
        y(1)=1/(A^2+B^2)*(B*C+A*sqrt((A^2+B^2)*R^2-C^2));
        y(2)=1/(A^2+B^2)*(B*C-A*sqrt((A^2+B^2)*R^2-C^2));
        x=-B/A*y+C/A;
    end

    Phi=atan2(y,x);					% Cylinder azimuth of intersection points
    psi=atan2((y-y_h),(x-x_h));		% Helix azimuth of intersection points
    phi=psi+h*pi/2;					% Track direction of intersection points

    % find the first intersection point
    n=[-2:2];
    n=n';
    deltaphi=phi-phi_s;
    deltaphin=[deltaphi(1)+n*2*pi,deltaphi(2)+n*2*pi];

    if h==1
        k1=find( deltaphin(:,1)>=0 & deltaphin(:,1)<2*pi );
        k2=find( deltaphin(:,2)>=0 & deltaphin(:,2)<2*pi );
    else
        k1=find( deltaphin(:,1)>-2*pi & deltaphin(:,1)<=0 );
        k2=find( deltaphin(:,2)>-2*pi & deltaphin(:,2)<=0 );
    end

    if abs(deltaphin(k1,1))<=abs(deltaphin(k2,2))
        deltaphi=deltaphin(k1,1);
        x=x(1);
        y=y(1);
		Phi=Phi(1);
    else
        deltaphi=deltaphin(k2,2);
        x=x(2);
        y=y(2);
        Phi=Phi(2);
	end

    z=z_s+h*Rh*cot(theta_s)*deltaphi;
    phi=phi_s+deltaphi;
    beta=phi-Phi;
    while beta>pi/2
        beta=beta-pi;
    end
    while beta<-pi/2
        beta=beta+pi;
    end
    if abs(sin(beta))>sinbmx
        ierr=2;
    end

    paramf=[Phi,z,theta_s,beta,kappa_s];
end