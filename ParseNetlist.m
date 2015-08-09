function [inputMap termMap instTermMap resValueMap capValueMap sizeG] = ParseNetlist(fileName)
  nlFile = fopen(fileName, 'r');
  line = fgets(nlFile);

  sizeG = 0;
  inputMap = containers.Map();
  termMap = containers.Map();
  instTermMap = containers.Map();
  resValueMap = containers.Map();
  capValueMap = containers.Map();
  while ischar(line)
    disp(line);
    if line(1) == '*'
      line = fgets(nlFile);
      continue;
    elseif line(1) == 'R'
      tokens = strsplit(line, ' ');
      instName = tokens(1);
      term1 = tokens(2);
      term2 = tokens(3);
      value = tokens(4);

      instName = instName{1};
      term1 = term1{1};
      term2 = term2{1};
      value = str2num(value{1});
      if strcmp(term1,'0') == 0 && isKey(termMap, term1) == 0
        termMap(term1) = length(termMap) + 1;
      end
      if strcmp(term2,'0') == 0 && isKey(termMap, term2) == 0
        termMap(term2) = length(termMap) + 1;
      end
      terms = {term1;term2};
      instTermMap(instName) = terms;
      resValueMap(instName) = value;
    elseif line(1) == 'C'
      term1 = tokens(1);
      term2 = tokens(2);
      value = tokens(3);
      term1 = term1{1};
      term2 = term2{1};
      value = str2num(value{1});
      if strcmp(term1,'0') == 0 && isKey(termMap, term1) == 0
        termMap(term1) = length(termMap) + 1;
      end
      if strcmp(term2,'0') == 0 && isKey(termMap, term2) == 0
        termMap(term2) = length(termMap) + 1;
      end
      terms = {term1;term2};
      instTermMap(instName) = terms;
      capValueMap(instName) = value;
    elseif line(1) == 'V'
      %Independent Voltage
      tokens = strsplit(line, ' ');
      cell = tokens(1);
      ivName = cell{1};
      cell = tokens(2);
      inputName1 = cell{1};
      cell = tokens(3);
      inputName2 = cell{1};

      cell = tokens(4);
      inputVoltage = cell{1};

      inputMap(inputName1) = str2num(inputVoltage);
      if strcmp(inputName2,'0') == 0 && isKey(termMap, inputName2) == 0
        termMap(inputName2) = length(termMap) + 1;
      end
    elseif line(1) == 'M'
      %MOS 
      tokens = strsplit(line, ' ');
      cell = tokens(1);
      instName = cell{1}
      cell = tokens(6);
      mosType = cell{1}
      mosParam = strjoin(tokens(7:length(tokens)), ' ')
      CreateInst(instName, mosType, mosParam);

      cell = tokens(2);
      d = cell{1}
      cell = tokens(3);
      g = cell{1}
      cell = tokens(4);
      s = cell{1}
      cell = tokens(5);
      b = cell{1}

      terms = {d;g;s;b};
      instTermMap(instName) = terms;
      for i = 1:length(terms)
        term = terms(i)
        term = term{1}
        if strcmp(term,'0') == 0 && isKey(termMap, term) == 0
          termMap(term) = length(termMap) + 1;
        end
      end
    end
    line = fgets(nlFile);
  end
  keys(termMap)
  keys(inputMap) 
  sizeG = length(termMap) + length(inputMap);
  nTerms = length(termMap);
  nInputs = length(inputMap);

  %%% initState = 0;
  %%% if X == 0 
  %%%   initState = 1;
  %%%   X = zeros(1,nTerms+nInputs)';
  %%% end
  %%% 
  %%% W = zeros(1,nTerms + nInputs)';
  %%% i = 1;
  %%% for input = keys(inputMap)
  %%%   input = input{1};
  %%%   W(nTerms + i) = inputMap(input);

  %%%   if initState == 1
  %%%     X(i) = inputMap(input);
  %%%   end

  %%%   i = i + 1;
  %%% end
  %%% W

  %%% G = [ zeros(nTerms,nTerms) eye(nTerms,nInputs) ;
  %%%       eye(nInputs,nTerms) zeros(nInputs, nInputs) ]

  %%% C = zeros(nTerms+nInputs, nTerms+nInputs)

  %%% instNameList = keys(instTermMap) 
  %%% for idx = 1:length(instNameList)
  %%%   instName = instNameList(idx);
  %%%   instName = instName{1};
  %%%   terms = instTermMap(instName)

  %%%   % Get Node Voltage
  %%%   V = [0;0;0;0];
  %%%   for i = 1:length(terms)
  %%%     termName = terms(i);
  %%%     termName = termName{1}
  %%%     if strcmp(termName, '0') == 0 
  %%%       index = termMap(termName);
  %%%       if isKey(inputMap, termName)
  %%%          V(i) = inputMap(termName);
  %%%       else
  %%%          if initState == 1
  %%%            V(i) = 0;
  %%%          else
  %%%            V(i) = X(index);
  %%%          end
  %%%       end
  %%%     end
  %%%   end
  %%%   [Gm, Cm, F, I, J, I2] = BSIMCMG(instName,V(1),V(2),V(3),V(4));
  %%%   %V
  %%%   X(nTerms+1:nTerms+nInputs) = zeros(1,nInputs)';
  %%%   
  %%%   rhs = zeros(1,nTerms+nInputs)';
  %%%   j = 1;
  %%%   for i = 1:length(terms)
  %%%     termName = terms(i);
  %%%     termName = termName{1}
  %%%     if strcmp(termName, '0') == 0 
  %%%       index = termMap(termName)
  %%%       G(index,index) = G(index,index) + -Gm(i,i);
  %%%       C(index,index) = C(index,index) + -Cm(i,i);
  %%%                       
  %%%       if isKey(termMap, termName)
  %%%          rhs(index) = rhs(index) + -F(i);
  %%%       end
  %%%       
  %%%       if isKey(inputMap, termName)
  %%%          X(nTerms+index) = X(nTerms+index) + -I(i);
  %%%       end
  %%%     end
  %%%   end
  %%%   %X(1:nTerms) = V(1:nTerms);

  %%%   rhs = rhs + W
  %%% end
  fclose(nlFile);
end

