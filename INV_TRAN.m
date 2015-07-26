function Y = INV_TRAN(Vdd, error, alpha)

  outfile = fopen('INV_TRAN.csv', 'w');

  Initialize();
  CreateInst('MN', 'nmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
  CreateInst('MP', 'pmos1', 'L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0');
    
  G = [ zeros(3,3) eye(3,2);
        eye(2,3) zeros(2,2)
  ];
       
  C = [ zeros(5,5) ];
  
  % Time step
  h = 1e-8
  % Time domain, from 0 to 2ms, time step is h 
  time = 0:h:1.5e-6
  % Store the BE result
  Y = zeros(5,length(time));
  
  %%%Backward Euler Method
  X0 = INV_DC(vin(time(1)), Vdd, error, alpha);
  Vout = X0(3);
  [PGm, PCm, PF, PI, PJ, PI2] = BSIMCMG('MP',Vdd,vin(time(1)),Vout,Vdd);
  [NGm, NCm, NF, NI, NJ, NI2] = BSIMCMG('MN',Vout,vin(time(1)),0,0);
  jp = -2.0/h * [PJ(1) + PJ(4); PJ(2) + NJ(2); PJ(3) + NJ(1); 0; 0];
  Xp = X0;
  Y(1:5,1) = Xp(1:5);
  fprintf(outfile, 'i,  t, Vdd, Vin, Vout, Ivdd, Ivin\n');
  fprintf(outfile, ['%d, %20g, %20g, %20g, %20g, %20g, %20g', '\n'], 1, time(1), Xp(1:5));
  fclose(outfile);
  
  for i = 2:length(time)
    outfile = fopen('INV_TRAN.csv', 'a');
    
    % Get DC operation point
    X = INV_DC(vin(time(i)), Vdd, error, alpha);
    Vout = X(3);  
    

    W = [0; 0; 0; Vdd; vin(time(i))];
    
    % For each element, get element MNA G, C
    [PGm, PCm, PF, PI, PJ, PI2] = BSIMCMG('MP',Vdd,vin(time(i)),Vout,Vdd);
    [NGm, NCm, NF, NI, NJ, NI2] = BSIMCMG('MN',Vout,vin(time(i)),0,0);
    
    % Update system MNA G, C
    G(1,1) = -PGm(1,1);
    G(2,2) = -PGm(2,2) + -NGm(2,2);
    G(3,3) = -PGm(3,3) + -NGm(1,1);
 
    C(1,1) = -PCm(1,1);
    C(2,2) = -PCm(2,2) + -NCm(2,2);
    C(3,3) = -PCm(3,3) + -NCm(1,1);

    M = ( G + 1/h * C );
    Pr = PI + PF - PGm*[Vdd;vin(time(i));Vout;Vdd];
    Nr = NI + NF - NGm*[Vout;vin(time(i));0;0];
    %r = I + F - Gm*[Vdd;vin(time(i));0;0];
    r = [Pr(1) + Pr(4); Pr(2) + Nr(2); Pr(3) + Nr(1); 0; 0];
    fn = [-PF(1) + -PF(4); -PF(2) + -NF(2); -PF(3) + -NF(1); 0; 0] + [r(1);r(2);r(3);0;0];
    jn = -1/h*[PJ(1) + PJ(4); PJ(2) + NJ(2); PJ(3) + NJ(1);0;0];

    % Update W
    W = [0; 0; 0; Vdd; vin(time(i))] + fn + jn - jp;
    
    % solve (G+1/h*C)Xn=1/h*C*Xp+W
    Xn = M\(1/h*C*Xp+W);
  
    % Output result
    Y(1:5,i) = Xn(1:5);
    
    fprintf(outfile, ['%d, %20g, %20g, %20g, %20g, %20g, %20g', '\n'], i, time(i), Xn(1:5));
    fclose(outfile);
    
    % Keep current state
    jp = jn;
    Xp = Xn;
  end
end




