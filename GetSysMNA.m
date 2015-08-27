function [G C X W rhs rhsF rhsJ] = GetSysMNA(X, inputMap, termMap, instTermMap, resValueMap, capValueMap, sizeG)
  %%% nlFile = fopen('input.sp', 'r');
  %%% line = fgets(nlFile);
  %%% sizeG = 0;
  %%% inputMap = containers.Map();
  %%% termMap = containers.Map();
  %%% instTermMap = containers.Map();
  %%% resValueMap = containers.Map();
  %%% capValueMap = containers.Map();
  %%% while ischar(line)
  %%%   disp(line);
  %%%   if line(1) == '*'
  %%%     line = fgets(nlFile);
  %%%     continue;
  %%%   elseif line(1) == 'R'
  %%%     tokens = strsplit(line, ' ');
  %%%     instName = tokens(1);
  %%%     term1 = tokens(2);
  %%%     term2 = tokens(3);
  %%%     value = tokens(4);

  %%%     instName = instName{1};
  %%%     term1 = term1{1};
  %%%     term2 = term2{1};
  %%%     value = str2num(value{1});
  %%%     if strcmp(term1,'0') == 0 && isKey(termMap, term1) == 0
  %%%       termMap(term1) = length(termMap) + 1;
  %%%     end
  %%%     if strcmp(term2,'0') == 0 && isKey(termMap, term2) == 0
  %%%       termMap(term2) = length(termMap) + 1;
  %%%     end
  %%%     terms = {term1;term2};
  %%%     instTermMap(instName) = terms;
  %%%     resValueMap(instName) = value;
  %%%   elseif line(1) == 'C'
  %%%     term1 = tokens(1);
  %%%     term2 = tokens(2);
  %%%     value = tokens(3);
  %%%     term1 = term1{1};
  %%%     term2 = term2{1};
  %%%     value = str2num(value{1});
  %%%     if strcmp(term1,'0') == 0 && isKey(termMap, term1) == 0
  %%%       termMap(term1) = length(termMap) + 1;
  %%%     end
  %%%     if strcmp(term2,'0') == 0 && isKey(termMap, term2) == 0
  %%%       termMap(term2) = length(termMap) + 1;
  %%%     end
  %%%     terms = {term1;term2};
  %%%     instTermMap(instName) = terms;
  %%%     capValueMap(instName) = value;
  %%%   elseif line(1) == 'V'
  %%%     %Independent Voltage
  %%%     tokens = strsplit(line, ' ');
  %%%     cell = tokens(1);
  %%%     ivName = cell{1};
  %%%     cell = tokens(2);
  %%%     inputName1 = cell{1};
  %%%     cell = tokens(3);
  %%%     inputName2 = cell{1};

  %%%     cell = tokens(4);
  %%%     inputVoltage = cell{1};

  %%%     inputMap(inputName1) = str2num(inputVoltage);
  %%%     if strcmp(inputName2,'0') == 0 && isKey(termMap, inputName2) == 0
  %%%       termMap(inputName2) = length(termMap) + 1;
  %%%     end
  %%%   elseif line(1) == 'M'
  %%%     %MOS 
  %%%     tokens = strsplit(line, ' ');
  %%%     cell = tokens(1);
  %%%     instName = cell{1}
  %%%     cell = tokens(6);
  %%%     mosType = cell{1}
  %%%     mosParam = strjoin(tokens(7:length(tokens)), ' ')
  %%%     CreateInst(instName, mosType, mosParam);

  %%%     cell = tokens(2);
  %%%     d = cell{1}
  %%%     cell = tokens(3);
  %%%     g = cell{1}
  %%%     cell = tokens(4);
  %%%     s = cell{1}
  %%%     cell = tokens(5);
  %%%     b = cell{1}

  %%%     terms = {d;g;s;b};
  %%%     instTermMap(instName) = terms;
  %%%     for i = 1:length(terms)
  %%%       term = terms(i)
  %%%       term = term{1}
  %%%       if strcmp(term,'0') == 0 && isKey(termMap, term) == 0
  %%%         termMap(term) = length(termMap) + 1;
  %%%       end
  %%%     end
  %%%   end
  %%%   line = fgets(nlFile);
  %%% end
  %%% keys(termMap)
  %%% keys(inputMap) 
  %%% sizeG = length(termMap) + length(inputMap);
  
  nTerms = length(termMap);
  nInputs = length(inputMap);

  %%initState = 0;
  %%if X == 0 
  %%  initState = 1;
  %%  X = zeros(1,nTerms+nInputs)';
  %%end
  
  W = zeros(1,nTerms + nInputs)';
  inputs = keys(inputMap);
  for i = 1:length(inputs)
    input = inputs(i);
    input = input{1};
    W(nTerms + i) = X(i);
  end
  W;

  G = [ zeros(nTerms,nTerms) eye(nTerms,nInputs) ;
        eye(nInputs,nTerms) zeros(nInputs, nInputs) ];

  C = zeros(nTerms+nInputs, nTerms+nInputs);

  %Update G for resistor
  disp('Update G for resistor');
  instRs = keys(resValueMap);
  for i = 1:length(instRs)
    instName = instRs(i);
    instName = instName{1};
    value = resValueMap(instName);
    terms = instTermMap(instName);
    term1 = terms(1);
    term1 = term1{1};
    term2 = terms(2);
    term2 = term2{1};
    if strcmp(term1, '0') == 1 && strcmp(term2, '0') == 0
      idxTerm2 = termMap(term2);
      G(idxTerm2, idxTerm2) = value;
    elseif strcmp(term1, '0') == 0 && strcmp(term2, '0') == 1
      idxTerm1 = termMap(term1);
      G(idxTerm1, idxTerm1) = value;
    elseif strcmp(term1, '0') == 0 && strcmp(term2, '0') == 0
      idxTerm2 = termMap(term2);
      idxTerm1 = termMap(term1);
      G(idxTerm1, idxTerm1) = value;
      G(idxTerm2, idxTerm2) = value;
      G(idxTerm1, idxTerm2) = -value;
      G(idxTerm2, idxTerm1) = -value;
    end
  end
  G

  %Update C for capacitor
  disp('Update C for resistor');
  instCs = keys(capValueMap);
  for i = 1:length(instCs)
    instName = instCs(i);
    instName = instName{1};
    value = capValueMap(instName);
    terms = instTermMap(instName);
    term1 = terms(1);
    term1 = term1{1};
    term2 = terms(2);
    term2 = term2{1};
    if strcmp(term1, '0') == 1 && strcmp(term2, '0') == 0
      idxTerm2 = termMap(term2);
      C(idxTerm2, idxTerm2) = value;
    elseif strcmp(term1, '0') == 0 && strcmp(term2, '0') == 1
      idxTerm1 = termMap(term1);
      C(idxTerm1, idxTerm1) = value;
    elseif strcmp(term1, '0') == 0 && strcmp(term2, '0') == 0
      idxTerm2 = termMap(term2);
      idxTerm1 = termMap(term1);
      C(idxTerm1, idxTerm1) = value;
      C(idxTerm2, idxTerm2) = value;
      C(idxTerm1, idxTerm2) = -value;
      C(idxTerm2, idxTerm1) = -value;
    end
  end
  C 
 
  X(nTerms+1:nTerms+nInputs) = zeros(1,nInputs)';
  rhs = zeros(1,nTerms+nInputs)';
  rhsF = zeros(1,nTerms+nInputs)';
  rhsJ = zeros(1,nTerms+nInputs)';
  rhs2 = zeros(1,nTerms+nInputs)';
    
  instNameList = keys(instTermMap);
  for idx = 1:length(instNameList)
    instName = instNameList(idx);
    instName = instName{1};
    terms = instTermMap(instName);

    % Get Node Voltage
    V = [0;0;0;0];
    for i = 1:length(terms)
      termName = terms(i);
      termName = termName{1};
      if strcmp(termName, '0') == 0 
        index = termMap(termName);
        %if isKey(inputMap, termName)
        %   %% V(i) = inputMap(termName);
        %    V(i) = X(index);
        %else
        %    V(i) = X(index);
        %end
        V(i) = X(index);
      end
    end

    if instName(1) == 'M'
      [Gm, Cm, F, I, J, I2] = BSIMCMG(instName,V(1),V(2),V(3),V(4))
      %V
      j = 1;
      for i = 1:length(terms)
        termName = terms(i);
        termName = termName{1};
        if strcmp(termName, '0') == 0 
          index = termMap(termName);
          G(index,index) = G(index,index) + -Gm(i,i);
          C(index,index) = C(index,index) + -Cm(i,i);
                          
          if isKey(termMap, termName)
             rhs(index) = rhs(index) + -F(i);
             rhsF(index) = rhsF(index) + F(i);
             rhsJ(index) = rhsJ(index) + J(i);
          end
          
          if isKey(inputMap, termName)
             idxInput = inputMap(termName);
             X(nTerms+idxInput) = X(nTerms+idxInput) + I(i);
             
             %X(nTerms+index) = X(nTerms+index) + I(i);
             
             %inputs = keys(inputMap);
             %for j = 1:length(inputs)
             %  input = inputs(j);
             %  input = input{1};
             %  if strcmp(input, termName) == 1
             %    X(nTerms + j) = X(nTerms + j) + I(i);
             %end
             
             rhs2(index) = rhs2(index) + -I(i);
          end
        end
      end
    elseif instName(1) == 'R'
      for i = 1:length(terms)
        term1 = terms(1);
        term2 = terms(2);
        term1 = term1{1};
        term2 = term2{1};
        
        vTerm1 = 0;
        vTerm2 = 0;
        if strcmp(term1, '0') == 0
          index = termMap(term1);
          vTerm1 = X(index);
        end
        if strcmp(term2, '0') == 0
          index = termMap(term2);
          vTerm2 = X(index);
        end

        Ir = (vTerm1-vTerm2)/resValueMap(instName); 
        
        if isKey(inputMap, term1)
          idxInput = inputMap(term1);
          X(nTerms+idxInput) = Ir;
        end
      end 

    end
    
    %X(1:nTerms) = V(1:nTerms);
  end
  rhs = rhs + W + rhs2;

  %%% fclose(nlFile);
end

