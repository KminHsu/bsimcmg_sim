function [G C X W rhs rhsF rhsJ loop isDiverge errors] = dc(Vin, Vdd, error, alpha, reset, nextOrderCnt)
  % chart 3-2: initialize the BSIMCMG, and read MOS model card.
  Initialize();

  errors = [];
  
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
           X(i) = Vin; 
         end
         if strcmp(termName, 'Vdd') == 1
           X(i) = Vdd; 
         end
      else
         disp('Outpout Node Index');
         disp(i);
         convNodeIdxs = [convNodeIdxs, i];
      end
    end
  end
  X0 = X;
  
  %[G C X rhs] = GetSysMNA(X, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG);
  %if Vin == 0 && Vdd == 1
  %  alpha = 1e5
  %  nextOrderCnt = 40;    
  %else
  %  nextOrderCnt = 50;
  %end
  alpha0 = alpha;
  
  isDiverge = 0;
  loop = 1;
  maxIter = 500;
  for loop = 1:maxIter
    [G C X W rhs rhsF rhsJ] = GetSysMNA(X, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG);
    f = G*X - rhs;

    %b = -f;
    %[L U P] = lu(G);
    %d = P*b;
    %y = L\d;
    %dX = U\y;

    %b = -f*alpha;
    %[L U P] = lu(G);
    %d = P*b;
    %y = L\d;
    %dX = U\y;
    %X = X + dX;
    %if mod(loop,nextOrderCnt) == 0
    %  alpha = alpha * 2
    %end
    %for j = 1:length(X)
    %  if X(j) >= reset
    %    X(j) = 1;
    %  end 
    %end

    %alp = -10:1:10;
    %alp = power(10,-alp);
    %minAlp = 0;
    %minF = 1e20;          
    %XX = X;
    %for k = 1:length(alp)
    %  nf = 0;
    %  [GG CC XX WW rhsrhs rhsFF rhsJJ] = GetSysMNA(XX + alp(k)*dX, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG);
    %  ff = GG*XX - rhsrhs;
    %  for j = 1:length(convNodeIdxs)
    %    l = convNodeIdxs(j);
    %    nf = nf + norm(ff(l));  
    %  end
    %  if nf < minF
    %    isCont = 0;  
    %    for j = 1:length(XX)
    %      if (XX(j) + alp(k)*dX(j)) >= reset
    %        isCont = 1;
    %      end 
    %    end
    %    if isCont == 1
    %      continue;
    %    end
    %    minF = nf;
    %    minAlp = alp(k);
    %  end
    %end
    %X = X + minAlp*dX;

    %test2 
    %b = -f;
    %[L U P] = lu(G);
    %d = P*b;
    %y = L\d;
    %dX = U\y;
    %isStop = 1;
    %for j = 1:length(convNodeIdxs)
    %  l = convNodeIdxs(j);
    %  alp = -2:1:10;
    %  alp = power(10,alp);
    %  minAlp = 1;
    %  minF = 0;          
    %  for m = 1:length(convNodeIdxs)
    %    ll = convNodeIdxs(m);
    %    minF = minF + norm(f(ll));  
    %  end
    %  for k = 1:length(alp)
    %    XX = X;
    %    XX(l) = XX(l) + alp(k)*dX(l)
    %    [GG CC XX WW rhsrhs rhsFF rhsJJ] = GetSysMNA(XX, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG);
    %    ff = GG*XX - rhsrhs;
    %    nf = 0;
    %    for m = 1:length(convNodeIdxs)
    %      ll = convNodeIdxs(m);
    %      nf = nf + norm(ff(ll));  
    %    end
    %    if nf < minF
    %      X(l) = X(l) + alp(k)*X(l);
    %    end
    %  end
    %end 
    %X
    
    %%Test1 
    %%b = -f;
    %%[L U P] = lu(G);
    %%d = P*b;
    %%y = L\d;
    %%dX = U\y;
    %%alp = -1:1:8;
    %%alp = power(10,alp);
    %%minF = 0;
    %%minAlp = 0;
    %%for j = 1:length(convNodeIdxs)
    %%  l = convNodeIdxs(j);
    %%  minF = minF + norm(f(l));  
    %%end
    %%for k = 1:length(alp)
    %%  nf = 0;
    %%  XX = X;
    %%  XX = XX + alp(k)*dX;
    %%  [GG CC XX WW rhsrhs rhsFF rhsJJ] = GetSysMNA(XX, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG);
    %%  ff = GG*XX - rhsrhs;
    %%  
    %%  for j = 1:length(convNodeIdxs)
    %%    l = convNodeIdxs(j);
    %%    nf = nf + norm(ff(l));  
    %%  end
    %%  
    %%  if nf < minF
    %%    minF = nf;
    %%    minAlp = alp(k);
    %%    break;
    %%  end
    %%end
    %%alpha = minAlp;
    %%X = X + alpha*dX;
   
    %Demo 
    b = -f*alpha;
    [L U P] = lu(G);
    d = P*b;
    y = L\d;
    dX = U\y;
    X = X + dX;

        
    f = G*X - rhs;
    allConv = 1;
    for j = 1:length(convNodeIdxs)
      %if abs(f(convNodeIdxs(j))) >= error
      if abs(dX(convNodeIdxs(j))) >= error
        allConv = 0;
      end
      errors = [errors abs(dX(convNodeIdxs(j)))];  
    end
    if allConv == 1
      loop
      break;
    end
    
    if mod(loop,nextOrderCnt) == 0
      alpha = alpha * 2
    end
    if reset > 0 
      for j = 1:length(convNodeIdxs)
        if abs(f(convNodeIdxs(j))) > reset
          nextOrderCnt = nextOrderCnt + nextOrderCnt;
          X(convNodeIdxs(j)) = X0(convNodeIdxs(j));
          alpha = alpha0;
        end
      end
    end
  end

  if loop == maxIter
    isDiverge = 1;    
  end
  %X = [Vdd; Vin; Vout; NI(1) + NI(4); PI(2) + NI(2)];
  if isDiverge == 1
    disp('Diverge');   
  end
end

