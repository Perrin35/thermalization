OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5705559) q[0];
sx q[0];
rz(-0.67222995) q[0];
sx q[0];
rz(1.4752522) q[0];
rz(-2.1885459) q[1];
sx q[1];
rz(-0.16521984) q[1];
sx q[1];
rz(-1.2686977) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27086285) q[0];
sx q[0];
rz(-2.6755973) q[0];
sx q[0];
rz(-0.52865882) q[0];
rz(-0.95999865) q[2];
sx q[2];
rz(-0.75163254) q[2];
sx q[2];
rz(3.1267191) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42434537) q[1];
sx q[1];
rz(-1.0043036) q[1];
sx q[1];
rz(-2.9866108) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3399063) q[3];
sx q[3];
rz(-0.82726631) q[3];
sx q[3];
rz(-2.8184067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98422009) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(2.3483707) q[2];
rz(-2.4685517) q[3];
sx q[3];
rz(-1.0854191) q[3];
sx q[3];
rz(1.8719155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9840045) q[0];
sx q[0];
rz(-0.80261153) q[0];
sx q[0];
rz(0.60930914) q[0];
rz(2.9318103) q[1];
sx q[1];
rz(-0.79133004) q[1];
sx q[1];
rz(-0.9598859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5989089) q[0];
sx q[0];
rz(-0.090714648) q[0];
sx q[0];
rz(-0.51341052) q[0];
rz(3.0716607) q[2];
sx q[2];
rz(-2.0581823) q[2];
sx q[2];
rz(2.5068381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7916225) q[1];
sx q[1];
rz(-2.4020139) q[1];
sx q[1];
rz(1.4776358) q[1];
x q[2];
rz(-1.0632564) q[3];
sx q[3];
rz(-1.1843269) q[3];
sx q[3];
rz(1.8563351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3557055) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(-2.6500224) q[2];
rz(-0.0056886557) q[3];
sx q[3];
rz(-2.0523043) q[3];
sx q[3];
rz(0.51386851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.046535) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(-2.3989578) q[0];
rz(0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(1.0661485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2258573) q[0];
sx q[0];
rz(-1.7122147) q[0];
sx q[0];
rz(-1.3456324) q[0];
rz(-pi) q[1];
rz(-1.891648) q[2];
sx q[2];
rz(-1.3898456) q[2];
sx q[2];
rz(2.2448283) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6934476) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(0.13685302) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1769165) q[3];
sx q[3];
rz(-2.6163273) q[3];
sx q[3];
rz(0.21533404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8891958) q[2];
sx q[2];
rz(-0.19813457) q[2];
sx q[2];
rz(0.55257094) q[2];
rz(-2.6342454) q[3];
sx q[3];
rz(-2.0523968) q[3];
sx q[3];
rz(0.56940091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4813389) q[0];
sx q[0];
rz(-0.57732552) q[0];
sx q[0];
rz(1.8950155) q[0];
rz(1.4056816) q[1];
sx q[1];
rz(-2.3697) q[1];
sx q[1];
rz(2.4235922) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3221949) q[0];
sx q[0];
rz(-0.85032636) q[0];
sx q[0];
rz(2.2706991) q[0];
rz(-pi) q[1];
rz(1.9393607) q[2];
sx q[2];
rz(-2.5466726) q[2];
sx q[2];
rz(-1.0234317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64252428) q[1];
sx q[1];
rz(-0.82083265) q[1];
sx q[1];
rz(1.7700511) q[1];
rz(2.4627732) q[3];
sx q[3];
rz(-2.9390237) q[3];
sx q[3];
rz(-2.5819419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(-0.70449746) q[2];
rz(-2.7590175) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(2.2883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.660897) q[0];
sx q[0];
rz(-3.1013885) q[0];
sx q[0];
rz(-2.8327292) q[0];
rz(2.1924696) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(2.5905051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2386467) q[0];
sx q[0];
rz(-1.4184877) q[0];
sx q[0];
rz(2.7480844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91294841) q[2];
sx q[2];
rz(-1.5163053) q[2];
sx q[2];
rz(-1.2555497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70610395) q[1];
sx q[1];
rz(-1.4248253) q[1];
sx q[1];
rz(1.7273497) q[1];
rz(0.14679175) q[3];
sx q[3];
rz(-0.78189584) q[3];
sx q[3];
rz(1.1679389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71517313) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(2.5185781) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(-0.4272517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039336786) q[0];
sx q[0];
rz(-0.66348851) q[0];
sx q[0];
rz(-0.46184194) q[0];
rz(-2.6448008) q[1];
sx q[1];
rz(-1.3235612) q[1];
sx q[1];
rz(1.7740446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5493126) q[0];
sx q[0];
rz(-0.89087109) q[0];
sx q[0];
rz(-0.10217765) q[0];
rz(-pi) q[1];
rz(0.15256792) q[2];
sx q[2];
rz(-1.0786453) q[2];
sx q[2];
rz(0.14482611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8277487) q[1];
sx q[1];
rz(-1.5428468) q[1];
sx q[1];
rz(-1.8716783) q[1];
rz(1.3472137) q[3];
sx q[3];
rz(-0.63399613) q[3];
sx q[3];
rz(-0.29350212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10609145) q[2];
sx q[2];
rz(-2.0095299) q[2];
sx q[2];
rz(2.3512225) q[2];
rz(-0.60449374) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-2.7214971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919493) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(-0.41437909) q[0];
rz(-0.62037933) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(1.5588123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65711439) q[0];
sx q[0];
rz(-3.0006391) q[0];
sx q[0];
rz(0.64224215) q[0];
rz(-pi) q[1];
rz(-1.5528658) q[2];
sx q[2];
rz(-1.6758783) q[2];
sx q[2];
rz(1.1370657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5276416) q[1];
sx q[1];
rz(-1.5613341) q[1];
sx q[1];
rz(-1.6285614) q[1];
x q[2];
rz(1.956432) q[3];
sx q[3];
rz(-1.7117071) q[3];
sx q[3];
rz(-2.8711365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9862426) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(0.17779329) q[2];
rz(-2.6627461) q[3];
sx q[3];
rz(-2.0279453) q[3];
sx q[3];
rz(-0.068304121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9947522) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(0.11257182) q[0];
rz(-0.62537891) q[1];
sx q[1];
rz(-1.1275147) q[1];
sx q[1];
rz(2.2523527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5035985) q[0];
sx q[0];
rz(-1.6498008) q[0];
sx q[0];
rz(-0.043009506) q[0];
rz(-0.28381376) q[2];
sx q[2];
rz(-0.2760016) q[2];
sx q[2];
rz(-2.2021879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.016805033) q[1];
sx q[1];
rz(-1.8062002) q[1];
sx q[1];
rz(1.0339869) q[1];
rz(-pi) q[2];
rz(-1.8440644) q[3];
sx q[3];
rz(-2.8972761) q[3];
sx q[3];
rz(-0.48297802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2832977) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(0.27321401) q[2];
rz(-2.9122399) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6943618) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(2.2253775) q[0];
rz(0.18237309) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.1590385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834476) q[0];
sx q[0];
rz(-0.30320692) q[0];
sx q[0];
rz(-0.68443735) q[0];
rz(-pi) q[1];
rz(-0.49075895) q[2];
sx q[2];
rz(-3.033394) q[2];
sx q[2];
rz(-2.8645017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7024544) q[1];
sx q[1];
rz(-1.7299486) q[1];
sx q[1];
rz(1.252878) q[1];
rz(-pi) q[2];
rz(2.4501514) q[3];
sx q[3];
rz(-0.85606282) q[3];
sx q[3];
rz(-2.5987529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3850022) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(-2.8263367) q[2];
rz(2.5473525) q[3];
sx q[3];
rz(-0.88163328) q[3];
sx q[3];
rz(0.63410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88157982) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(-2.9823533) q[0];
rz(-2.0533994) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(-0.27096567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20734678) q[0];
sx q[0];
rz(-1.1109612) q[0];
sx q[0];
rz(0.28069541) q[0];
rz(-pi) q[1];
rz(3.0130544) q[2];
sx q[2];
rz(-2.4646799) q[2];
sx q[2];
rz(2.235379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9393255) q[1];
sx q[1];
rz(-1.9069873) q[1];
sx q[1];
rz(2.4656117) q[1];
x q[2];
rz(-0.48153523) q[3];
sx q[3];
rz(-1.2709444) q[3];
sx q[3];
rz(-1.7026342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3403885) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(0.32269746) q[2];
rz(3.0835551) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(2.8431852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67644453) q[0];
sx q[0];
rz(-1.6315176) q[0];
sx q[0];
rz(2.3254707) q[0];
rz(-2.8836518) q[1];
sx q[1];
rz(-2.0057269) q[1];
sx q[1];
rz(-1.534091) q[1];
rz(-0.28889523) q[2];
sx q[2];
rz(-1.8452273) q[2];
sx q[2];
rz(0.065963521) q[2];
rz(1.1318351) q[3];
sx q[3];
rz(-0.96228941) q[3];
sx q[3];
rz(-1.8214772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
