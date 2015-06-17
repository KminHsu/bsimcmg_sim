% h = 0.1
% X0 = [2 ;0 ;0 ;0 ;0]
% G = [ 1 -1 0 0 1 ; -1 1 0 1 0; 0 0 1 -1 0; 0 1 -1 0 0; 1 0 0 0 0]
% C = [ 0 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -2 0; 0 0 0 0 0 ]
% W = [0 ;0 ;0 ;0 ;1]

% G = [ 1 -1 1 ; -1 1 0; 1 0 0 ]
% C = [ 0 0 0; 0 1 0; 0 0 0 ]
% W = [0 ;0 ;1]

G = [ 1e-3 -1e-3 1 ; -1e-3 1e-3 0; 1 0 0 ] + [ 0 0 0; 0 1e-3 0; 0 0 0]
C = [ 0 0 0; 0 0.001e-6 0; 0 0 0 ]
h = 1e-8
x = 0:h:2e-6
BE_y = zeros(length(x),1);
TR_y = zeros(length(x),1);

%%%BE
W = [0; 0; vin(1)];
X0 = inv(G)*W
X = X0;
M = ( G + 1/h * C );
for i = 2:length(x)
  W = [0 ; 0; vin(x(i))];
  X = inv(M)*(1/h*C*X+W);
  BE_y(i) = X(2);
end

%%%TR
W = [0; 0; vin(1)];
X0 = inv(G)*W;
X = X0;
M1 = ( G + 2/h*C );
M2 = -1*( G - 2/h*C );
for i = 2:length(x)
  W = [0 ; 0; vin(x(i))];
  X = inv(M1)*(M2*X + W + W);
  TR_y(i) = X(2);
end

%figure
%plot(x, BE_y-TR_y)
figure
plot(x, BE_y)
%figure
%plot(x, TR_y, '*')
%plot(x, BE_y, '*', x, TR_y, 'O')

% %%%BE
% X = X0;
% M = ( G + 1/h * C );
% for i = 0:100
%   X = inv(M)*(1/h*C*X+W);
% end
% BE_X = X
% 
%  
% %%%TR
% X = X0;
% M1 = ( G + 2/h*C );
% M2 = -1*( G - 2/h*C );
% for i = 0:100
%   X = inv(M1)*(M2*X + W + W);
% end
% TR_X = X
% 
% %%%N-R & BE
% X = X0;
% M = ( G + 1/h * C );
% for i = 0:5
%   F = M*X - (1/h*C*X+W);
%   J = M-1/h*C;
%   %dX = -1*inv(J)*F;
%   dX = -1*J\F;
%   X = X + dX;
% end
% NR_BE_X = X
% 
% %%%N-R & BE
% X = X0;
% M1 = ( G + 2/h*C );
% M2 = -1*( G - 2/h*C );
% for i = 0:5
%   F = M1*X - (M2*X + W + W);
%   J = M1 - M2;
%   %dX = -1*inv(J)*F;
%   dX = -1*J\F;
%   X = X + dX;
% end
% NR_TR_X = X

