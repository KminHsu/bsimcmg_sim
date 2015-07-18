% nMOS
function X = MOS_DC(Vd,Vg,Vs,Vb)
  Initialize();
  CreateInst('MN', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
 
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
  [NGm, NQm, NF, NI, NJ, NI2] = BSIMCMG('MN',Vd,Vg,Vs,Vb);
  % update the MNA matrix, G.
  for row = 1:4
    G(row,1:4) = -NGm(row,1:4);
  end
  r = NI + NF - NGm*[Vd;Vg;Vs;Vb];
  % solve GX=W
  X = G\ ([-NF(1); -NF(2); -NF(3); -NF(4); 0; 0; 0 ; 0] + W + [r(1);r(2);r(3);r(4);0;0;0;0]);
end
