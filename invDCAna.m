% INV DC Analysis
Vdd = 0:0.1:0;
Vin = 0:0.1:1.0;
MId = zeros(length(Vin), length(Vdd));
for i = 1:length(Vdd)
  for j = 1:length(Vin)
    X = INV_DC(Vin(j),Vdd(i));
    MId(i,j) = X(3);
  end
end
%plot(MId);

outfile = fopen('INV_DC.csv', 'w')
fprintf(outfile, ['Vin,' repmat('%20g,',1,length(Vin)) '\n'], Vin)
fprintf(outfile, ['Vdd,' repmat('%20g,',1,length(Vdd)) '\n'], Vdd)
for row = 1:length(Vin)
  fprintf(outfile, [repmat('%20g,',1,length(Vdd)) '\n'], MId(row,:));
end
fclose(outfile);
