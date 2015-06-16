% h = 0.1
% X0 = [2 ;0 ;0 ;0 ;0]
% G = [ 1 -1 0 0 1 ; -1 1 0 1 0; 0 0 1 -1 0; 0 1 -1 0 0; 1 0 0 0 0]
% C = [ 0 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -2 0; 0 0 0 0 0 ]
% W = [0 ;0 ;0 ;0 ;1]

% G = [ 1 -1 1 ; -1 1 0; 1 0 0 ]
% C = [ 0 0 0; 0 1 0; 0 0 0 ]
% W = [0 ;0 ;1]

G = [ 1 -1 1 ; -1 1 0; 1 0 0 ]
C = [ 0 0 0; 0 1 0; 0 0 0 ]
W = [0 ;0 ;1]

% h = 0.001
h = 5e-6
len = 20000
x = linspace(0.0001, 0.1, len);
BE_y = linspace(0.0001, 0.1, len);
TR_y = linspace(0.0001, 0.1, len);

%%%BE
X0 = inv(G)*W
X = X0;
M = ( G + 1/h * C );
for i = 1:length(x)
  if i ==  floor(len/2)
    W(3) = 0
  end
  X = inv(M)*(1/h*C*X+W);
  BE_y(i) = X(2);
end
BE_X = X;

%%%TR
W = [0 ;0 ;1];
X0 = inv(G)*W
X = X0;
M1 = ( G + 2/h*C );
M2 = -1*( G - 2/h*C );
for i = 1:length(x)
  if i ==  floor(len/2)
    W(3) = 0
  end
  X = inv(M1)*(M2*X + W + W);
  TR_y(i) = X(2);
end
TR_X = X;

figure
plot(x, BE_y-TR_y)
figure
plot(x, BE_y)
figure
plot(x, TR_y, '*')

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

