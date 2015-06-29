% h = 0.1
% X0 = [2 ;0 ;0 ;0 ;0]
% G = [ 1 -1 0 0 1 ; -1 1 0 1 0; 0 0 1 -1 0; 0 1 -1 0 0; 1 0 0 0 0]
% C = [ 0 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -2 0; 0 0 0 0 0 ]
% W = [0 ;0 ;0 ;0 ;1]

% G = [ 1 -1 1 ; -1 1 0; 1 0 0 ]
% C = [ 0 0 0; 0 1 0; 0 0 0 ]
% W = [0 ;0 ;1]


Initialize();
CreateInst('M1', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');

G = [ 0  0  0  0  1  0 ; 
      0  0  0  0  0  1 ;
      0  0  0  0  0  0 ;
      0  0  0  0  0  0 ;
      1  0  0  0  0  0 ;
      0  1  0  0  0  0 ];
     
C = [ 0  0  0  0  0  0 ; 
      0  0  0  0  0  0 ;
      0  0  0  0  0  0 ;
      0  0  0  0  0  0 ;
      0  0  0  0  0  0 ;
      0  0  0  0  0  0 ];

[Gm, Qm, F, I] = BSIMCMG('M1',0,0,0,0)

for i = 1:4
  for j = 1:4
    G(i + (j-1)*6) = Gm(i + (j-1)*4);
    C(i + (j-1)*6) = Qm(i + (j-1)*4);
  end
end

% Time step
h = 1e-8
% Time domain, from 0 to 2ms, time step is h
x = 0:h:2e-6
% Store the BE result
BE_y = zeros(length(x),1);
% store the TR result
TR_y = zeros(length(x),1);

%%%Backward Euler Method
W = [0; 0; 0; 0; 0.3; 1.0];
[Gm, Qm, F, I] = BSIMCMG('M1',0.3, 1.0, 0.0, 0.0)
X0 = [ 0.3; 1.0; 0; 0; I(1); I(2)]
X = X0;
for i = 2:length(x)
  [Gm, Qm, F, I] = BSIMCMG('M1',X(1), X(2), X(3), X(4))
  for i = 1:4
    for j = 1:4
      G(i + (j-1)*6) = Gm(i + (j-1)*4);
      C(i + (j-1)*6) = Qm(i + (j-1)*4);
    end
  end

  M = ( G + 1/h * C );
  W = [0; 0; 0; 0; X(1); X(2)];
  % In here, I don't use Newton method. I just use inverse function to solve this equation
  X = inv(M)*(1/h*C*X+W);
  % Store Vout, that is Vc
  BE_y(i) = X(5);
end

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
