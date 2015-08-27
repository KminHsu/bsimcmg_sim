function tran(Vdd, error, alpha, reset, nextOrderCnt, outFileName)

  outfile = fopen(outFileName, 'w');

  Initialize();
  
  % Time step
  h = 1e-8
  % Time domain, from 0 to 2ms, time step is h 
  time = 0:h:1.5e-6
  
  [inputMap termMap instTermMap resValueMap capValueMap sizeG] = ParseNetlist('input.sp');

  nTerms = length(termMap);
  nInputs = length(inputMap);
  
  convNodeIdxs = [];
  
  X = zeros(1,nTerms+nInputs)';
  X0 = zeros(1,nTerms+nInputs)';
  terms = keys(termMap);
  for i = 1:length(terms)
    termName = terms(i);
    termName = termName{1}
    if strcmp(termName, '0') == 0 
      index = termMap(termName);
      if isKey(inputMap, termName)
         X(i) = inputMap(termName);
         if strcmp(termName, 'Vin') == 1
           X(i) = vin(time(1)); 
         end
         if strcmp(termName, 'Vdd') == 1
           X(i) = Vdd; 
         end
      else
         %disp('Outpout Node Index');
         %disp(i);
         convNodeIdxs = [convNodeIdxs, i];
      end
    end
  end
  X0 = X;
  

  
  %%%Backward Euler Method
  [G C X0 W rhs rhsF rhsJ loop isDiverge] = dc(vin(time(1)), Vdd, error, alpha, reset, nextOrderCnt);
  
  jp = -2.0/h * rhsJ;
  Xp = X0;
    
  fprintf(outfile, ['i,  t,' repmat('X%d,',1,length(Xp))  '\n'], 1:1:length(Xp));
  fprintf(outfile, ['%d, %20g,' repmat('%20g,',1,length(Xp)) '%20g, %20g', '\n'], 1, time(1), Xp, loop, isDiverge);
  fclose(outfile);
  
  for i = 2:length(time)
    outfile = fopen(outFileName, 'a');
    
    % Get DC operation point
    [G C X W rhs rhsF rhsJ loop isDiverge] = dc(vin(time(i)), Vdd, error, alpha, reset,nextOrderCnt);

    M = ( G + 1/h * C );
    %Pr = PI + PF - PGm*[Vdd;vin(time(i));Vout;Vdd];
    %Nr = NI + NF - NGm*[Vout;vin(time(i));0;0];
    %r = [Pr(1) + Pr(4); Pr(2) + Nr(2); Pr(3) + Nr(1); 0; 0];
    %fn = [-PF(1) + -PF(4); -PF(2) + -NF(2); -PF(3) + -NF(1); 0; 0] + [r(1);r(2);r(3);0;0];
    
    fn = -rhsF; %[-PF(1) + -PF(4); -PF(2) + -NF(2); -PF(3) + -NF(1); 0; 0];
    jn = -1/h*rhsJ; %[PJ(1) + PJ(4); PJ(2) + NJ(2); PJ(3) + NJ(1);0;0];

    % Update W1
    W1 = W + fn + jn - jp;
    
    % solve (G+1/h*C)Xn=1/h*C*Xp+W
    Xn = M\(1/h*C*Xp+W1);
      
    fprintf(outfile, ['%d, %20g,' repmat('%20g,',1,length(Xn)) '%20g, %20g', '\n'], i, time(i), Xn, loop, isDiverge);

    fclose(outfile);
    
    % Keep current state
    jp = jn;
    Xp = Xn;
  end
end




