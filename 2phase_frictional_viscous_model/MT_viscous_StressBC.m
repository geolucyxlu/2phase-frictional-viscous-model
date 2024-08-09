function [E_MT,M_MT,Em,Sm,e,s,vis,vis0] = MT_viscous_StressBC(Sigma,a,q,eta,eta0,Nc,Nm,n,phi)
%MORI_TANAKA_MULTIPHASE_STRESSBC Summary of this function goes here
%   Detailed explanation goes here
% 08/July/2021
%  matrix using Mori-Tanaka homogenization appraoch (nonlinear,various
%  shapes inclusions, isotropic matrix)
%  ! stress boundary condition !
%  Input Parameters:-------------------------------------------------------
% Sigma:Imposed macroscale stress tensor, a 3-by-3 matrix.
% n:	The number of inclusions.
% a:	The initial three semi-axes of clasts (a1>a2>a3; a3=1), a 3-by-3 matrix.
% q:	The initial orientations of clasts as transformation matrix, 
%       a 3-by-3-by-n matrix. 
% eta:	Clasts' effective viscosities defined at the initial state, 
%       such as the imposed macroscopic strain-rate state, a 1-by-n vector.
% eta0: Matrix viscosity defined at the imposed satate.
% Nc:	Clasts' stress exponents, a 1-by-n vector.
% Nm:	Matrix stress exponent.
% phi:  Clast volume fraction.
%% ------------------------------------------------------------------------ 
%  generate 4th-order identity tensors    
       [Jd, ~, ~, ~] = FourIdentity();
       
       % matrix stiffness
       Cm = 2*eta0*Jd;
       Mm = FourTensorInv(Cm);
       %  strain rate invariants at which the viscosities of RDEs are defiend  
       REF  =Inva(Sigma)*ones(1,n);  
       REF0 =Inva(Sigma);
       Sm = Sigma;
       Ce = zeros(3,3,3,3,n);
       Me = Ce;
      for j=1:n
          Ce(:,:,:,:,j) = 2*eta(j)*Jd;
          Me(:,:,:,:,j) = 1/(2*eta(j))*Jd;
      end
    
      B_dil = Ce;
      mB    = Ce;
      s = zeros(3,3,n); 
      e = s;
      %% Mori-Tanaka
     for kk=1:1000
         
         for r=1:n
              SS=SnPI_vis(a(:,r));
              %  transform the Eshelby tensor(SS, PI) to the global coordinate
              S=Transform(SS,q(:,:,r)');
              
             for jj=1:1000
                 %global coordinate
                 t1 = Contract(S,Contract(Mm,Ce(:,:,:,:,r))); %secant compliances
                 t2 = Jd-S+Nm*t1;
                 t3 = FourTensorInv(t2);
                 t4 = Jd+(Nm-1)*S;                                                  
                 A_dil= Contract(t3,t4);     %A in Eq.23 Jiang,2014
                 t5 = Contract(Ce(:,:,:,:,r),A_dil);
                 B_dil(:,:,:,:,r)=Contract(t5,Mm);
                 
                 % partition strain rate from matrix strain rate
                 s(:,:,r) = Multiply(B_dil(:,:,:,:,r),Sm); 
                 %     the second invariant of the strain rate tensor
                 sI = (0.5*contract1(s(:,:,r) ,s(:,:,r) ))^0.5;
                 zz       = (sI/REF(r))^(1-Nc(r))*eta(r) ;
%                  % check viscosity 
%                  alpha    = abs(zz/eta(r) - 1);
                 % check stress
                 alpha    = abs(sI/REF(r)-1);
                 eta(r)   = zz;
                 REF(r)   = sI;
                 Ce(:,:,:,:,r)  = 2*eta(r)*Jd; 
                 Me(:,:,:,:,r)  = 1/(2*eta(r))*Jd;
              
                 if alpha < 0.001 
                     break
                 end
                 if jj==1000
                        XX=['warning: clast ',num2str(r),' not converged!!'];
                        disp(XX);
                 end
             end
             e(:,:,r) = Multiply(Me(:,:,:,:,r) ,s(:,:,r));
             mB(:,:,:,:,r)=Contract(Me(:,:,:,:,r),B_dil(:,:,:,:,r));
         end
         BB = (1-phi)*Jd + phi*mean(B_dil,5);
         Bm_MT = FourTensorInv(BB);
         Sm = Multiply(Bm_MT,Sigma);
         % the second invariant of the strain rate tensor 
         SmI = (0.5*contract1(Sm,Sm))^0.5;
         % update viscosity 
         yy  = (SmI/REF0)^(1-Nm)*eta0;
%          % check viscosiy
%          beta    = abs(yy/eta0 - 1);
         % check stress
         beta    = abs(SmI/REF0 - 1);
         eta0    = yy;
         REF0    = SmI;
         Cm = 2*eta0*Jd; 
         Mm = FourTensorInv(Cm);
         
         if  beta < 0.001
             break
         end
         if kk==1000
                YY='warning: Matrix is not converged!!';
                disp(YY);
         end   
     end
        
        vis=eta;
        vis0=eta0;
        mB1 = (1-phi)*Mm + phi*mean(mB,5);
        M_MT = Contract(mB1,Bm_MT);
        E_MT = Multiply(M_MT,Sigma);
        Em   = Multiply(Mm,Sm);

end

