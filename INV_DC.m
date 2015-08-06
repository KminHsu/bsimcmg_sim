function [X loop isDiverge] = INV_DC(Vin, Vdd, error, alpha)
  % chart 3-2: initialize the BSIMCMG, and read MOS model card.
  Initialize();
  CreateInst('MN', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  CreateInst('MP', 'pmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  % chart 3-3: initialize matrices and variables.
  G = [ zeros(3,3) eye(3,2);
        eye(2,3) zeros(2,2)
  ];
       
  C = [ zeros(5,5) ];
  
  % Time step
  h = 1e-8
  % Time domain, from 0 to 2ms, time step is h
  x = 0:h:1.5e-6
  
  %%%Backward Euler Method
  %X = [0;0;0;0;0];
  %if Vin == 0.0 && Vdd == 0.0
  % return
  %end
 
  %initial guss
  Vout = 0;
  %if Vdd == 1.0 && Vin == 1.0 
  %  Vout =  3.05807731e-08;
  %elseif Vdd == 1.0 && Vin == -1.0
  %  Vout = 1.0;
  %elseif Vdd == 1.0 && Vin == 0.0
  %  Vout =  9.99987736e-01;
  %else  
  %  Vout = 10e-8;
  %end
  alpha0 = alpha;
  isDiverge = 0;
  loop = 1;
  nextOrderCnt = 5;
  maxIter = 2000;
  for loop = 1:maxIter
    W = [0; 0; 0; Vdd; Vin];
    % chart 3-4-2: Calculate MNA matrices, Gm, Qm, F, I, J, I2 by Newton’s iteration methods.
    [PGm, PQm, PF, PI, PJ, PI2] = BSIMCMG('MP',Vdd,Vin,Vout,Vdd);
    [NGm, NQm, NF, NI, NJ, NI2] = BSIMCMG('MN',Vout,Vin,0,0);
    % chart 3-4-1: update the MNA matrix, G.
    G(1,1) = -PGm(1,1) + -PGm(4,4);
    G(2,2) = -PGm(2,2) + -NGm(2,2);
    G(3,3) = -PGm(3,3) + -NGm(1,1);
  
    Pr = PI + PF - PGm*[Vdd;Vin;Vout;Vdd];
    Nr = NI + NF - NGm*[Vout;Vin;0;0];
    rhs = ([-PF(1) + -PF(4); -PF(2) + -NF(2); -PF(3) + -NF(1); 0; 0] + W + -[ NI(1) + NI(4); PI(2) + NI(2); 0; 0; 0] );
    f = G*[Vdd; Vin; Vout; PI(3) + NI(1); PI(2) + NI(2) ] - rhs;
    % chart 3-5: ||Xk+1 - Xk|| < error?
    % f(3)  is Vout
    if abs(f(3)) < error
      loop
      break
    end
    % chart 3-6: update the damped Newton’s method to solve FX = GX + W.
    %dVout = -alpha*100*f(3)/G(3,3);
    
    %dX = G\-f*alpha;
    %dVout = dX(3);
    
    b = -f*alpha;
    [L U P] = lu(G);
    d = P*b;
    y = L\d;
    dX = U\y;
    dVout = dX(3);

    % chart 3-7: Xn+1(t)=Xn(t). 
    Vout = Vout + dVout;
    
    if mod(loop,nextOrderCnt) == 0
      alpha = alpha * 2
    end
    
    if Vout > 5.0 || Vout < -5.0
      %nextOrderCnt = nextOrderCnt * 2;
      %nextOrderCnt = nextOrderCnt * nextOrderCnt;
      nextOrderCnt = nextOrderCnt + nextOrderCnt;
      Vout = 0;
      alpha = alpha0;
    end
  end
  if loop == maxIter
    isDiverge = 1;    
  end
  X = [Vdd; Vin; Vout; NI(1) + NI(4); PI(2) + NI(2)];
  if isDiverge == 1
    disp('Diverge');   
  end
end

