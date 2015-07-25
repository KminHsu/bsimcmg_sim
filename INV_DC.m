function X = INV_DC(Vin, Vdd, error, alpha)
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
  % Store the BE result
  BE_y1 = zeros(length(x),1);
  BE_y2 = zeros(length(x),1);
  BE_y3 = zeros(length(x),1);
  BE_y4 = zeros(length(x),1);
  BE_y5 = zeros(length(x),1);
  BE_y6 = zeros(length(x),1);
  
  %%%Backward Euler Method
  %Vdd = 0.1;
  %Vout = 3.73947213e-02;
  
  if Vin == 0.0 && Vdd == 0.0
   X = [0;0;0;0;0]
   return
  end
 
  %initial guss 
  Vout = 10e-8;
  for loop = 1:20000
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
    dVout = -alpha*100*f(3)/G(3,3);
    % chart 3-7: Xn+1(t)=Xn(t). 
    Vout = Vout + dVout;
    
    if mod(loop,200) == 0
      alpha = alpha * 10
    end
    
    if Vout > 2.0
      disp('Diverge');
      loop
      break
    end
    
  end
  loop
  X = [Vdd; Vin; Vout; NI(1) + NI(4); PI(2) + NI(2)];
end

% fn = -[F(1);F(2);F(3);F(NI(1) + NI(4); PI(2) + NI(2)4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
% jp = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];

% BE_y1(1) = Xp(1);
% BE_y2(1) = Xp(2);
% BE_y3(1) = Xp(3);
% BE_y4(1) = Xp(4);
% BE_y5(1) = Xp(5);
% BE_y6(1) = Xp(6);
% BE_y7(1) = Xp(7);
% BE_y8(1) = Xp(8);

% for i = 2:length(x)
%   [Gm, Qm, F, I, J, I2] = BSIMCMG('M1',-0.3, -vin(x(i)),0,0);
%   for k = 1:4
%     for j = 1:4
%       G(k + (j-1)*8) = -Gm(k + (j-1)*4);
%       C(k + (j-1)*8) = -Qm(k + (j-1)*4);
%     end
%   end
%   M = ( G + 1/h * C );
%   r = I + F - Gm*[-0.3;-vin(x(i));0;0];
%   fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
%   jn = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];
% 
%   W = [0; 0; 0; 0; -0.3; -vin(x(i)); 0; 0] + fn + jn - jp;
%   %W = [0; 0; 0.3; vin(x(i))] + fn + jp;
%   %W = [0; 0; 0.3; vin(x(i))] + fn;
%   Xn = M\(1/h*C*Xp+W);
% 
%   BE_y1(i) = Xn(1);
%   BE_y2(i) = Xn(2);
%   BE_y3(i) = Xn(3);
%   BE_y4(i) = Xn(4);
%   BE_y5(i) = Xn(5);
%   BE_y6(i) = Xn(6);
%   BE_y7(i) = Xn(7);
%   BE_y8(i) = Xn(8);
%   
%   jp = jn;
%   Xp = Xn;
% 
% end
% plot(x, BE_y4);

% %%%Trapezoido Method
% Wp = [0; 0; 0.3; vin(x(1))];
% [Gm, Qm, F, I, J, I2] = BSIMCMG('M1',0.3,vin(x(1)),0,0);
% for k = 1:2
%   for j = 1:2
%     G(k + (j-1)*4) = -Gm(k + (j-1)*4);
%     C(k + (j-1)*4) = -Qm(k + (j-1)*4);
%   end
% end
% rp = I + F - Gm*[0.3;vin(x(1));0;0;];
% fp = [-F(1);-F(2);0;0] + [rp(1); rp(2); 0; 0];
% jp = [-2.0/h*J(1);-2.0/h*J(2);0;0];
% X0 = G\ ( Wp + fp );
% Xp = X0;
% 
% BE_y1(1) = Xp(1);
% BE_y2(1) = Xp(2);
% BE_y3(1) = Xp(3);
% BE_y4(1) = Xp(4);
% BE_I(1) = I(2);
% 
% for i = 2:length(x)
%   [Gm, Qm, F, I, J, I2] = BSIMCMG('M1',0.3, vin(x(i)), 0, 0);
%   for k = 1:2
%     for j = 1:2
%       G(k + (j-1)*4) = -Gm(k + (j-1)*4);
%       C(k + (j-1)*4) = -Qm(k + (j-1)*4);
%     end
%   end
% 
%   M = ( G + 2.0/h*C );
% 
%   Wn = [0; 0; 0.3; vin(x(i))];
%   rn = I + F - Gm*[0.3;vin(x(i));0;0;];
%   fn = [-F(1);-F(2);0;0] + [rn(1); rn(2); 0; 0];
%   jn = [-2.0/h*J(1);-2.0/h*J(2);0;0];
% 
%   Xn = M \ ( -(G-2.0/h*C)*Xp  + Wn + Wp + fn + jn + fp - jp );
% 
%   BE_y1(i) = Xn(1);
%   BE_y2(i) = Xn(2);
%   BE_y3(i) = Xn(3);
%   BE_y4(i) = Xn(4);
%   BE_I(i) = I(2);
% 
%   Wp = Wn;
%   fp = fn;
%   jp = jn;
%   Xp = Xn;
%   rp = rn;
% end
% plot(x, BE_y4);



%-----------------------------------------------------------------------
% % There are two resistors
% G = [ 1e-3 -1e-3 1 ; -1e-3 1e-3 0; 1 0 0 ] + [ 0 0 0; 0 1e-3 0; 0 0 0]
% % Only one capacitor
% C = [ 0 0 0; 0 0.001e-6 0; 0 0 0 ]

% % Time step
% h = 1e-8
% % Time domain, from 0 to 2ms, time step is h
% x = 0:h:2e-6
% % Store the BE result
% BE_y = zeros(length(x),1);
% % store the TR result
% TR_y = zeros(length(x),1);

% %%%Backward Euler Method
% W = [0; 0; vin(1)];
% X0 = inv(G)*W
% X = X0;
% M = ( G + 1/h * C );
% for i = 2:length(x)
%   W = [0 ; 0; vin(x(i))];
%   % In here, I don't use Newton method. I just use inverse function to solve this equation
%   X = inv(M)*(1/h*C*X+W);
%   % Store Vout, that is Vc
%   BE_y(i) = X(2);
% end
% 
% %%%Trapzoidal Method
% W = [0; 0; vin(1)];
% X0 = inv(G)*W;
% X = X0;
% M1 = ( G + 2/h*C );
% M2 = -1*( G - 2/h*C );
% for i = 2:length(x)
%   W = [0 ; 0; vin(x(i))];
%   % In here, I don't use Newton method. I just use inverse function to solve this equation
%   X = inv(M1)*(M2*X + W + W);
%   % Store Vout, that is Vc
%   TR_y(i) = X(2);
% end
% 
% %figure
% %plot(x, BE_y-TR_y)
% figure
% plot(x, BE_y)
% %figure
% %plot(x, TR_y, '*')
% %plot(x, BE_y, '*', x, TR_y, 'O')



