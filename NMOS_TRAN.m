% nMOS
function [Y, time] = NMOS_TRAN(Vdd) 
  % chart 2-2: initialize the BSIMCMG, and read MOS model card.
  Initialize();
  CreateInst('MN', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  % chart 2-3: initialize matrices and variables.
  % circuit MNA matrix
  G = [ 0  0  0  0  1  0  0  0; 
        0  0  0  0  0  1  0  0;
        0  0  0  0  0  0  1  0;
        0  0  0  0  0  0  0  1;
        1  0  0  0  0  0  0  0;
        0  1  0  0  0  0  0  0;
        0  0  1  0  0  0  0  0;
        0  0  0  1  0  0  0  0;
  ];
       
  C = [ 0  0  0  0  0  0  0  0;   
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0
  ];
  
  % Time step
  h = 1e-8
  % Time domain, from 0 to 2ms, time step is h
  time = 0:h:1.5e-6
  % Store the BE result
  Y = zeros(8, length(time));
  
  %%%Backward Euler Method
  W = [0; 0; 0; 0; Vdd; vin(time(1)); 0 ; 0];
  [Gm, Qm, F, I, J, I2] = BSIMCMG('MN',Vdd,vin(time(1)),0,0);
  
  G(1:4,1:4) = -Gm(1:4,1:4);
  C(1:4,1:4) = -Qm(1:4,1:4);
  
  r = I + F - Gm*[Vdd;vin(time(1));0;0];
  % chart 2-4: using the Backward Euler method to solve the X0 (t = 0, DC analysis).
  % solve GX0=W.
  X0 = G\ ([-F(1); -F(2); -F(3); -F(4); 0; 0; 0 ; 0] + W + [r(1);r(2);r(3);r(4);0;0;0;0]);
  Xp = X0;
  fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
  jp = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];
  
  Y(1:8,1) = Xp(1:8);
  % chart 2-5: increasing t.
  for i = 2:length(time)
    % chart 2-6-2: calculate MNA matrices, Gm, Qm, F, I, J, I2 by Newton’s iteration methods.
    [Gm, Qm, F, I, J, I2] = BSIMCMG('MN',Vdd, vin(time(i)),0,0);
    % chart 2-6-1: update the MNA matrices, G, C, and W.
    G(1:4,1:4) = -Gm(1:4,1:4);
    C(1:4,1:4) = -Qm(1:4,1:4);

    M = ( G + 1/h * C );
    r = I + F - Gm*[Vdd;vin(time(i));0;0];
    fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
    jn = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];
    % chart 2-7: prepare W.
    W = [0; 0; 0; 0; Vdd; vin(time(i)); 0; 0] + fn + jn - jp;
    
    % chart 2-8: solve (G+1/h*C)Xn=1/h*C*Xp+W.
    Xn = M\(1/h*C*Xp+W);
  
    Y(1:8,i) = Xn(1:8);
    % chart 2-9: Xn(t+1) = Xn(t).
    jp = jn;
    Xp = Xn;
  end
end
