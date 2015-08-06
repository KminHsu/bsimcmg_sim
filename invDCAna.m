% INV DC Analysis
Vdd = 0:0.1:1.0;
Vin = 0:0.1:1.0;
MVout = zeros(length(Vin), length(Vdd));

outfile = fopen('INV_DC.log', 'w')
error = 1e-18;
alpha = 1e4;
for i = 1:length(Vdd)
  for j = 1:length(Vin)
    outfile = fopen('INV_DC.log', 'a');
    [X loop isDiverge] = INV_DC(Vin(j),Vdd(i),error,alpha);
    MVout(i,j) = X(3);
    fprintf(outfile, ['%g %g %g %g %g', '\n'], Vin(j), Vdd(i), X(3), loop, isDiverge);
    fclose(outfile);
  end
end

% error = 1e-18;
% alpha = 1;
% i = 1;
% for j = 1:length(Vin)
%   X = INV_DC(Vin(j),Vdd(i), error, alpha);
%   MVout(i,j) = X(3);
% end

outfile = fopen('INV_DC.csv', 'w')
fprintf(outfile, ['Vin,' repmat('%20g,',1,length(Vin)) '\n'], Vin)
fprintf(outfile, ['Vdd,' repmat('%20g,',1,length(Vdd)) '\n'], Vdd)
for row = 1:length(Vin)
  fprintf(outfile, [repmat('%20g,',1,length(Vdd)) '\n'], MVout(row,:));
end
fclose(outfile);
