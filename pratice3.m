


% C = [ 0 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -2 0; 0 0 0 0 0 ]
% W = [0 ;0 ;0 ;0 ;1]

% G = [ 1 -1 1 ; -1 1 0; 1 0 0 ]
% C = [ 0 0 0; 0 1 0; 0 0 0 ]
% W = [0 ;0 ;1]


Initialize();
CreateInst('M1', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');

G = [ 0  0  0  0  1  0  0  0; 
      0  0  0  0  0  1  0  0;
      0  0  0  0  0  0  1  0;
      0  0  0  0  0  0  0  1;
      1  0  0  0  0  0  0  0;
      0  1  0  0  0  0  0  0;
      0  0  1  0  0  0  0  0;
      0  0  0  1  0  0  0  0;
];
     
C = [ 0  0  0  0  0  0  0  0;   
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0
];

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
BE_y7 = zeros(length(x),1);
BE_y8 = zeros(length(x),1);

%%%Backward Euler Method
W = [0; 0; 0; 0; 0.3; vin(x(1)); 0 ; 0];
[Gm, Qm, F, I, J, I2] = BSIMCMG('M1',0.3,vin(x(1)),0,0);
for k = 1:4
  for j = 1:4
    G(k + (j-1)*8) = -Gm(k + (j-1)*4);
    C(k + (j-1)*8) = -Qm(k + (j-1)*4);
  end
end
r = I + F - Gm*[0.3;vin(x(1));0;0];
X0 = G\ ([-F(1); -F(2); -F(3); -F(4); 0; 0; 0 ; 0] + W + [r(1);r(2);r(3);r(4);0;0;0;0]);
Xp = X0;
fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
jp = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];

BE_y1(1) = Xp(1);
BE_y2(1) = Xp(2);
BE_y3(1) = Xp(3);
BE_y4(1) = Xp(4);
BE_y5(1) = Xp(5);
BE_y6(1) = Xp(6);
BE_y7(1) = Xp(7);
BE_y8(1) = Xp(8);

for i = 2:length(x)
  [Gm, Qm, F, I, J, I2] = BSIMCMG('M1',0.3, vin(x(i)),0,0);
  for k = 1:4
    for j = 1:4
      G(k + (j-1)*8) = -Gm(k + (j-1)*4);
      C(k + (j-1)*8) = -Qm(k + (j-1)*4);
    end
  end
  M = ( G + 1/h * C );
  r = I + F - Gm*[0.3;vin(x(i));0;0];
  fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
  jn = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];

  W = [0; 0; 0; 0; 0.3; vin(x(i)); 0; 0] + fn + jn - jp;
  %W = [0; 0; 0.3; vin(x(i))] + fn + jp;
  %W = [0; 0; 0.3; vin(x(i))] + fn;
  Xn = M\(1/h*C*Xp+W);

  BE_y1(i) = Xn(1);
  BE_y2(i) = Xn(2);
  BE_y3(i) = Xn(3);
  BE_y4(i) = Xn(4);
  BE_y5(i) = Xn(5);
  BE_y6(i) = Xn(6);
  BE_y7(i) = Xn(7);
  BE_y8(i) = Xn(8);
  
  jp = jn;
  Xp = Xn;

end
plot(x, BE_y4);

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

