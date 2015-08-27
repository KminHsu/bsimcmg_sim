% INV DC Analysis
Vdd = 0:0.1:1.0;
Vin = 0:0.1:1.0;
MVout = zeros(length(Vin), length(Vdd));

outfile = fopen('INV_DC.log', 'w')
fclose(outfile);
error = 1e-12;
alpha = 1e6;
reset = 10;
nextOrderCnt = 60;
for i = 1:length(Vdd)
  for j = 1:length(Vin)
    outfile = fopen('INV_DC.log', 'a');
    %[G C X W rhs loop isDiverge] = INV_DC_2(Vin(j),Vdd(i),error,alpha);
    [G C X W rhs rhsF rhsJ loop isDiverge errors] = dc(Vin(j), Vdd(i), error, alpha, reset, nextOrderCnt)
    MVout(i,j) = X(3);
    fprintf(outfile, ['%g %g %g %g %g', '\n'], Vin(j), Vdd(i), X(3), loop, isDiverge);
    fprintf(outfile, ['errors->' repmat('%20g,',1,length(errors)) '\n'], errors);
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
