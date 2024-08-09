function [S_MT,C_MT,Em,Sm,e,s,vis,vis0] = MT_viscous(E,a,q,eta,eta0,Nc,Nm,n,phi)
% Mori_Tanaka.m
%   07/07/2021     
%  description: To get bulk rehology for composite material with a connected
%  matrix using Mori-Tanaka homogenization appraoch (nonlinear,various
%  shapes inclusions, isotropic matrix)
%  ! strain rate boundary condition !
%  Input Parameters:-------------------------------------------------------
% E:	Imposed macroscale strain rate tensor, a 3-by-3 matrix.
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
   REF  = (0.5* contract1(E,E))^0.5*ones(1,n);  
   REF0 =(0.5* contract1(E,E))^0.5;
   Em = E;
   Ce = zeros(3,3,3,3,n);
  for j=1:n
      Ce(:,:,:,:,j) = 2*eta(j)*Jd;
  end

  A_dil = Ce;
  %Ar_MT = Ce;
  cA    = Ce;
  e = zeros(3,3,n); 
  s = e;
  %%Mori-Tanaka
 for kk=1:1000
     
     parfor r=1:n

          SS=SnPI_vis(a(:,r));
          %  transform the Eshelby tensor(SS, PI) to the global coordinate
          S=Transform(SS,q(:,:,r)');
          
         for jj=1:1000
             %global coordinate
             t1 = Contract(S,Contract(Mm,Ce(:,:,:,:,r))); %secant compliances
             t2 = Jd-S+Nm*t1;
             t3 = FourTensorInv(t2);
             t4 = Jd+(Nm-1)*S;                                                  
             A_dil(:,:,:,:,r) = Contract(t3,t4);     %A in Eq.23 Jiang,2014
             % partition strain rate from matrix strain rate
             e(:,:,r) = Multiply(A_dil(:,:,:,:,r),Em); 
             %     the second invariant of the strain rate tensor
             epsilonI = (0.5*contract1(e(:,:,r) ,e(:,:,r) ))^0.5;
             zz       = (epsilonI/REF(r))^((1-Nc(r))/Nc(r))*eta(r) ;
             alpha    = abs(zz/eta(r) - 1);
             eta(r)   = zz;
             REF(r)   = epsilonI;
             Ce(:,:,:,:,r)  = 2*eta(r)*Jd; 
             s(:,:,r) = Multiply(Ce(:,:,:,:,r),e(:,:,r));
             if alpha < 0.001 
                 break
             end
             if jj==1000
                    XX=['warning: clast ',num2str(r),' not converged!!'];
                    disp(XX);
             end
         end
         cA(:,:,:,:,r)=Contract(Ce(:,:,:,:,r),A_dil(:,:,:,:,r));
     end
     AA = (1-phi)*Jd + phi*mean(A_dil,5);
     Am_MT = FourTensorInv(AA);
     Em= Multiply(Am_MT,E);
     % the second invariant of the strain rate tensor 
     epsilonI0 = (0.5*contract1(Em,Em))^0.5;
     % update viscosity 
     yy       = (epsilonI0/REF0)^((1-Nm)/Nm)*eta0;
     % check strain rate
     beta    = abs(epsilonI0/REF0 - 1);
     eta0     = yy;
     REF0      = epsilonI0;
     Cm  = 2*eta0*Jd; 
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
    cA1 = (1-phi)*Cm + phi*mean(cA,5);
    C_MT = Contract(cA1,Am_MT);
    S_MT = Multiply(C_MT,E);
    Sm = Multiply(Cm,Em);
end

   

