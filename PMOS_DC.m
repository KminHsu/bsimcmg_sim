% pMOS
function X = MOS_DC(Vd,Vg,Vs,Vb)
  % chart 1-2: initialize the BSIMCMG, and read MOS model card.
  Initialize();
  CreateInst('MP', 'pmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  % chart 1-3: initialize matrices and variables.
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
  
  W = [0; 0; 0; 0; Vd;Vg;Vs;Vb];
  [PGm, PQm, PF, PI, PJ, PI2] = BSIMCMG('MP',Vd,Vg,Vs,Vb);
  % update the MNA matrix, G
  for row = 1:4
    G(row,1:4) = -PGm(row,1:4);
  end
  r = PI + PF - PGm*[Vd;Vg;Vs;Vb];
  % chart 1-4: solve GX=W.
  X = G\ ([-PF(1); -PF(2); -PF(3); -PF(4); 0; 0; 0 ; 0] + W + [r(1);r(2);r(3);r(4);0;0;0;0]);
end
