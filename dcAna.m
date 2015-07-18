% NMOS DC Analysis
Vd = 0:0.1:1.0;
Vg = 0:0.1:1.5;
MId = zeros(length(Vd), length(Vg));
for i = 1:length(Vd)
  for j = 1:length(Vg)
    X = NMOS_DC(Vd(i),Vg(j),0,0);
    MId(i,j) = X(5);
  end
end
%plot(MId);

outfile = fopen('MOS_DC.csv', 'w')
fprintf(outfile, ['Vg,' repmat('%20g,',1,length(Vg)) '\n'], Vg)
fprintf(outfile, ['Vd,' repmat('%20g,',1,length(Vd)) '\n'], Vd)
for row = 1:length(Vd)
  fprintf(outfile, [repmat('%20g,',1,length(Vg)) '\n'], MId(row,:));
end
%fclose(outfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMOS DC Analysis
MId = zeros(length(Vd), length(Vg));
for i = 1:length(Vd)
  for j = 1:length(Vg)
    X = PMOS_DC(-Vd(i),-Vg(j),0,-Vd(i));
    MId(i,j) = X(5);
  end
end
%plot(MId);

%outfile = fopen('MOS_DC.csv', 'w')
fprintf(outfile, ['Vg,' repmat('%20g,',1,length(Vg)) '\n'], -Vg)
fprintf(outfile, ['Vd,' repmat('%20g,',1,length(Vd)) '\n'], -Vd)
for row = 1:length(Vd)
  fprintf(outfile, [repmat('%20g,',1,length(Vg)) '\n'], MId(row,:));
end
fclose(outfile);




