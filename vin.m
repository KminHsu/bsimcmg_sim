function y = vin(t)
  if t == 0.0
    y = 0.0;
  elseif t >= 1e-7 && t <= 5e-7
    y = 5.0;
  else
    y = 0.0;
  end
end
