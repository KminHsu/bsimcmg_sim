** RC Network, high pass
** tran(1.0, 1e-18, 1, 10, 50, 'rc_h.csv')
*Vvin Vin 0 0.3
*C1 Vin n1 100e-12
*R1 n1 0 1000

** RC Network, low pass
** tran(1.0, 1e-18, 1, 10, 50, 'rc_l.csv')
Vvin Vin 0 0.3
R1 Vin n1 1000
C1 n1 0 100e-12


** 1 inverter
** tran(1.0, 1e-12, 1e6, 10, 40, 'inv2.csv')
*Vvdd Vdd 0 1.0
*Vvin Vin 0 0.3
*MP1 Vdd Vin Vout1 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN1 Vout1 Vin 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0

** 12 inverter path
** tran(1.0, 1e-18, 1e6, 10, 50, 'inv12.csv')
*Vvdd Vdd 0 1.0
*Vvin Vin 0 0.3
*MP1 Vdd Vin Vout1 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN1 Vout1 Vin 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP2 Vdd Vout1 Vout2 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN2 Vout2 Vout1 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP3 Vdd Vout2 Vout3 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN3 Vout3 Vout2 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP4 Vdd Vout3 Vout4 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN4 Vout4 Vout3 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP5 Vdd Vout4 Vout5 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN5 Vout5 Vout4 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP6 Vdd Vout5 Vout6 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN6 Vout6 Vout5 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP7 Vdd Vout6 Vout7 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN7 Vout7 Vout6 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP8 Vdd Vout7 Vout8 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN8 Vout8 Vout7 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP9 Vdd Vout8 Vout9 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN9 Vout9 Vout8 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP10 Vdd Vout9 Vout10 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN10 Vout10 Vout9 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP11 Vdd Vout10 Vout11 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN11 Vout11 Vout10 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MP12 Vdd Vout11 Vout12 Vdd pmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0
*MN12 Vout12 Vout11 0 0 nmos1 L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0

** MOSAmp
** tran(1.0, 1e-12, 1, 10, 1, 'mosamp.csv')
*Vvdd Vdd 0 1.0
*Vvin Vin 0 0.3
*R1 Vdd Vout 1000
*C1 Vout 0 10e-12
*MN1 Vout Vin 0 0 nmos1 L=6e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0


