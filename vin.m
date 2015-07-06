% function y = vin(t)
%   if t == 0.0
%     y = 0.0;
%   elseif t > 0 && t <= 5e-7
%     %y = 5.0;
%     y = 1.0;
%   else
%     y = 0.0;
%   end
% end

function y = vin(t)
  h = 1e-8;
  if t == 0.0
    y = 0.0;
  elseif t >= h && t <= 1e-7
    x1 = h;
    x2 = 1e-7;
    y1 = 0.0;
    y2 = 1.0;
    y = (t-x1) * (y2-y1)/(x2-x1) + y1;
  elseif t > 1e-7 && t <= 5e-7
    y = 1.0;
  elseif t >= 5e-7 && t <= 6e-7
    x1 = 5e-7;
    x2 = 6e-7;
    y1 = 1.0;
    y2 = 0.0;
    y = (t-x1) * (y2-y1)/(x2-x1) + y1;
  else
    y = 0.0;
  end
end

