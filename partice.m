h = 0.1
X0 = [2 ;0 ;0 ;0 ;0]
G = [ 1 -1 0 0 1 ; -1 1 0 1 0; 0 0 1 -1 0; 0 1 -1 0 0; 1 0 0 0 0]
C = [ 0 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -2 0; 0 0 0 0 0 ]
W = [0 ;0 ;0 ;0 ;1]

h = 0.0001
x = linspace(0.0001, 0.1, 100000)
BE_y = linspace(0.0001, 0.1, 100000)
TR_y = linspace(0.0001, 0.1, 100000)

%%%BE
X = X0;
M = ( G + 1/h * C );
for i = 1:length(x)
  X = inv(M)*(1/h*C*X+W);
  BE_y(i) = X(3);
end
BE_X = X;

%%%TR
X = X0;
M1 = ( G + 2/h*C );
M2 = -1*( G - 2/h*C );
for i = 1:length(x)
  X = inv(M1)*(M2*X + W + W);
  TR_y(i) = X(3);
end
TR_X = X;

figure
%plot(x, BE_y, '*', x, TR_y, 'O')
plot(x, BE_y-TR_y)

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

