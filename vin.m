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
  if t == 0.0
    y = 0.0;
  elseif t >= h && t <= 1e-7
    %y = 5.0;
  elseif t>= h && t <= 1e-7
    y = t * 1.0/(1e-7-h) + H;
  else
    y = 0.0;
  end
end

