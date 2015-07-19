% pMOS
function [Y, time] = PMOS_TRAN(Vdd) 
  Initialize();
  CreateInst('MP', 'pmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  
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
  W = [0; 0; 0; 0; Vdd; -vin(time(1)); 0 ; Vdd];
  [Gm, Qm, F, I, J, I2] = BSIMCMG('MP',Vdd,-vin(time(1)), 0 ,Vdd);
  
  G(1:4,1:4) = -Gm(1:4,1:4);
  C(1:4,1:4) = -Qm(1:4,1:4);

  r = I + F - Gm*[Vdd;-vin(time(1));0;Vdd];
  X0 = G\ ([-F(1); -F(2); -F(3); -F(4); 0; 0; 0 ; 0] + W + [r(1);r(2);r(3);r(4);0;0;0;0]);
  Xp = X0;
  fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
  jp = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];
  
  Y(1:8,1) = Xp(1:8);
  
  for i = 2:length(time)
    [Gm, Qm, F, I, J, I2] = BSIMCMG('MP',Vdd, -vin(time(i)),0,Vdd);

    G(1:4,1:4) = -Gm(1:4,1:4);
    C(1:4,1:4) = -Qm(1:4,1:4);

    M = ( G + 1/h * C );
    r = I + F - Gm*[Vdd;-vin(time(i));0;Vdd];
    fn = -[F(1);F(2);F(3);F(4);0;0;0;0] + [r(1);r(2);r(3);r(4);0;0;0;0];
    jn = -1/h*[J(1);J(2);J(3);J(4);0;0;0;0];

    W = [0; 0; 0; 0; Vdd; -vin(time(i)); 0; Vdd] + fn + jn - jp;
    
    % solve (G+1/h*C)Xn=1/h*C*Xp+W
    Xn = M\(1/h*C*Xp+W);
  
    Y(1:8,i) = Xn(1:8);
    
    jp = jn;
    Xp = Xn;
  end
end
