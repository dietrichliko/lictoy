function [paramf,dermat,alrphi,ierr]=prop5mod(radi,parami,idir,radf,zmin,zmax,sinbmx,iopt)
%                                                                     
%    aim :                                                             
%    -----                                                             
%    extrapolate a helix defined by the initial parameters parami      
%    up to a given cylinder, and compute if requested the derivatives  
%    of the final parameters w.r.t. the initial ones                   
%                                                                      
%    the computation uses double precision on intermediate variables   
%    if the variation of phi angle is less than dphimn (.0001 in this  
%    version) the computation is done at first order in 1/r in order   
%    to avoid rounding errors, especially in the derivatives           
%                                                                      
%    input  :  radi          : initial radius
%    input  :  parami(1-5)   : initial parameters                      
%                              (Phi,z,theta,beta,1/r)                
%                              with beta = Phi-phi                     
%                                       geometrical sign)              
%              idir    :  if  1 : positive extrapolation only          
%                         if -1 : negative         "                   
%                         if  0 : extrapolation on both sides          
%              radf          : radius of the cylinder                  
%              zmin          : lower z limit of the cylinder           
%              zmax          : upper z limit of the cylinder           
%              sinbmx        : maximum allowed for |sin(beta)| at the  
%                              intersection                            
%              iopt          : 0 if derivatives not requested          
%                              1 if derivatives requested              
%                                                                      
%    output :  ierr          : 0 if intersection found                 
%                              1 if no intersection with the cylinder  
%                              2 if sinbmx exceeded                    
%                              3 if intersection outside of limits     
%              paramf(1-5)   : final parameters                        
%              dermat(5,5)   : deriv. of final w.r.t. initial param.   
%                              der(1) = d(phi)/d(theta)                
%                              der(2) = d(phi)/d(beta)                 
%                              der(3) = d(phi)/d(1/r)                  
%                              der(4) = d(z)/d(theta)                  
%                              der(5) = d(z)/d(beta)                   
%                              der(6) = d(z)/d(1/r)                    
%                              der(7) = d(beta)/d(beta)                
%                              der(8) = d(beta)/d(1/r)                 
%              alrphi        : length (in r-phi projection) from start 
%                              to extrapolation, with a sign (positive 
%                              if the extrapolation is towards the     
%                              direction defined by theta,phi)         
%                                                                      
%    author  :  p. billoir                                             
%                                                                      
%    first version : 26-01-88                                          
%                                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
paramf=zeros(5,1);
dermat=zeros(5,5);
alrphi=0;
twopi=pi*2;
dphimn=1e-4;
%
ierr=0;
%
dermat=eye(5);
%
phii=parami(1);
beta=parami(4);
cosb=cos(beta);
sinb=sin(beta);
cotth=1./tan(parami(3));
rtrk=1./parami(5);
%
%  center and squared radius of the projected helix, in a rotated frame
%  (x-axis through the starting point)
xc=radi-rtrk*sinb;
yc=rtrk*cosb;
rc2=xc^2+yc^2;
%
%  resolution of a second order equation
rrr=(radf^2-rtrk^2-rc2)/(2.*rtrk);
delt=rc2-rrr^2;
%
%   intersection exists if delt > 0
%
if delt<=0
   ierr=1;
else
   delt=sqrt(delt);
   %
   %   choose intersection on the same side as the starting point
   %   w.r.t. the plane containing the z axis and the axis of the helix
   %
   sinf=(xc*rrr+yc*delt)/rc2;
   cosf=(xc*delt-yc*rrr)/rc2;
   xf=xc+rtrk*sinf;
   yf=yc-rtrk*cosf;
   sinbf=(sinf*xf-cosf*yf)/radf;
   %
   %   exit if beta too large at the intersection
   if abs(sinbf)>sinbmx
      ierr=2;
   else
      alff=atan2(sinf,cosf);
      dphi=alff-beta;
      alrphi=rtrk*dphi;
      %
      %   select positive or negative extrapolations, or both
      %
      if alrphi*idir<0
         ierr=1;
         return
      end
      %
      %   switch to approximate expressions if the variation of phi
      %   is less than dphimn
      %
      %   "exact" expressions ---------------------------------------
      %
      if abs(dphi)>=dphimn
         %
         zf=parami(2)+cotth*rtrk*dphi;
         %
         %   exit if outside of limits in z,theta,phi
         tth=pi/2.-atan(zf/radf);
         phif=atan2(yf,xf);
         pph=phif;
         if (pph<0) 
            pph=pph+twopi;
         end
         if (zf<zmin | zf>zmax)
            ierr=3;
         end
         %
         %   final parameters
         paramf(1,1)=mod(phii+phif+twopi,twopi);
         paramf(2,1)=zf;
         paramf(3,1)=parami(3);
         paramf(4,1)=alff-phif;
         paramf(5,1)=parami(5);
         %
         %   computation of derivatives ---------
         %
         if iopt==1
             cosbf=sqrt(1.-sinbf^2);
             %
             %  ccpsi = rc*cos(cap.psi) ; scpsi = rc*sin(cap.psi)     (initial point)
             %  ccpsf = rc*cos(cap.psi) ; scpsf = rc*sin(cap.psi)     (final point)
             ccpsi=radi-rtrk*sinb;
             scpsi=rtrk*cosb;
             ccpsf=radf-rtrk*sinbf;
             scpsf=rtrk*cosbf;
             %
             %  cpsii = sgn*rc*cos(psi) ; spsii = sgn*rc*sin(psi)     (initial point)
             %  cpsif = sgn*rc*cos(psi) ; spsif = sgn*rc*sin(psi)     (final point)
             cpsii=rtrk-radi*sinb;
             spsii=-radi*cosb;
             cpsif=rtrk-radf*sinbf;
             spsif=-radf*cosbf;
             %
             sdphi=sin(dphi);
             cdphi=cos(dphi);
             %
             fact=-rtrk/spsif;
             dermat(1,4)=sdphi*fact;
             dermat(1,5)=fact*rtrk*(1.-cdphi);
             dermat(2,3)=-rtrk*dphi*(1.+cotth^2);
             dermat(2,4)=rtrk*cotth*(radf*ccpsf*spsii/spsif-radi*ccpsi)/rc2;
             dermat(2,5)=rtrk^2*cotth*(-dphi+sinbf/cosbf-...
                 (radi*scpsi+radf*ccpsf*cpsii/spsif)/rc2);
             dermat(4,4)=spsii/spsif;
             dermat(4,5)=rtrk*(cpsif-cpsii)/spsif;
         end

         %
         %   approximation at first order in 1/r --------------------------
         %
      else
         dr2=radf^2-radi^2;
         rcosb=radi*cosb;
         aa=1.-radi*sinb/rtrk;
         delt=rcosb^2+aa*dr2;
         %   exit if no solution
         if delt<=0
            ierr=1;
         else
            rdphi=(sqrt(delt)-rcosb)/aa;
            dphi=rdphi/rtrk;
            dcosf=-sinb-.5*cosb*dphi;
            cosf=cosb+dcosf*dphi;
            yf=-rdphi*dcosf;
            dsinf=cosb-.5*sinb*dphi;
            sinf=sinb+dsinf*dphi;
            xf=radi+rdphi*dsinf;
            sinbf=(sinf*xf-cosf*yf)/radf;
            zf=parami(2)+cotth*rdphi;
            %
            %
            %   exit if outside of limits in z,theta,phi
            tth=pi/2.-atan(zf/radf);
            phif=atan2(yf,xf);
            pph=phif;
            if pph<0 
               pph=pph+twopi;
            end
            if (zf<zmin | zf>zmax) 
               ierr=3;
               %keyboard
            end
            phif=atan2(yf,xf);
            paramf(1,1)=mod(parami(1)+phif+twopi,twopi);
            paramf(2,1)=zf;
            paramf(3,1)=parami(3);
            paramf(4,1)=parami(4)+dphi-phif;
            paramf(5,1)=parami(5);
            %
            %   computation of derivatives --------------
            %
            if iopt==1
                cosbf=sqrt(1.-sinbf^2);
                sphif=yf/radf;
                %
                dermat(1,4)=rdphi/(radf*cosbf);
                dermat(1,5)=.5*rdphi*dermat(1,4);
                dermat(2,3)=-rdphi*(1.+cotth^2);
                dermat(2,4)=radi*cotth*sphif/cosbf;
                dermat(2,5)=.5*rdphi*dermat(2,4);
                dermat(4,4)=(radi*cosb)/(radf*cosbf);
                dermat(4,5)=.5*rdphi*(1.+dermat(4,4));
            end
         end
      end
   end
end
%
return
