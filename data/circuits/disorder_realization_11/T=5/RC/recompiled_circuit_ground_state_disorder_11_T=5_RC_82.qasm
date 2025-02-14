OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.017555822) q[0];
sx q[0];
rz(-2.9889844) q[0];
sx q[0];
rz(0.18327644) q[0];
rz(-2.2493989) q[1];
sx q[1];
rz(-2.0018938) q[1];
sx q[1];
rz(-2.5633864) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59435049) q[0];
sx q[0];
rz(-0.092012398) q[0];
sx q[0];
rz(3.0334183) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(0.13730857) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44408724) q[1];
sx q[1];
rz(-1.1050535) q[1];
sx q[1];
rz(2.5986141) q[1];
x q[2];
rz(1.6434604) q[3];
sx q[3];
rz(-2.6748219) q[3];
sx q[3];
rz(0.49589402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9649428) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(-3.0435666) q[2];
rz(0.81884223) q[3];
sx q[3];
rz(-3.0374073) q[3];
sx q[3];
rz(-2.5067743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1408511) q[0];
sx q[0];
rz(-2.8241557) q[0];
sx q[0];
rz(0.66824085) q[0];
rz(0.60403281) q[1];
sx q[1];
rz(-0.030345358) q[1];
sx q[1];
rz(-2.277453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3383822) q[0];
sx q[0];
rz(-0.59640455) q[0];
sx q[0];
rz(-2.9687702) q[0];
rz(0.42704849) q[2];
sx q[2];
rz(-2.2042254) q[2];
sx q[2];
rz(2.065393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64795463) q[1];
sx q[1];
rz(-1.3820547) q[1];
sx q[1];
rz(2.8961839) q[1];
x q[2];
rz(-0.20012466) q[3];
sx q[3];
rz(-1.311655) q[3];
sx q[3];
rz(2.5166496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16022564) q[2];
sx q[2];
rz(-1.1796494) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(-1.9085599) q[3];
sx q[3];
rz(-2.8630856) q[3];
sx q[3];
rz(-1.2109141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8104372) q[0];
sx q[0];
rz(-1.5758608) q[0];
sx q[0];
rz(-2.1203777) q[0];
rz(1.637623) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(-2.5655897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5344328) q[0];
sx q[0];
rz(-2.7078621) q[0];
sx q[0];
rz(-2.8612616) q[0];
x q[1];
rz(-1.0025315) q[2];
sx q[2];
rz(-0.45323786) q[2];
sx q[2];
rz(1.1627878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29167029) q[1];
sx q[1];
rz(-1.909316) q[1];
sx q[1];
rz(2.4531219) q[1];
rz(-1.3542451) q[3];
sx q[3];
rz(-1.3050021) q[3];
sx q[3];
rz(2.1590421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7716498) q[2];
sx q[2];
rz(-2.6910431) q[2];
sx q[2];
rz(3.0453299) q[2];
rz(0.43109584) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(-2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189321) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(-0.2952964) q[0];
rz(0.88116208) q[1];
sx q[1];
rz(-0.22717871) q[1];
sx q[1];
rz(-0.49240246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8618362) q[0];
sx q[0];
rz(-1.4040213) q[0];
sx q[0];
rz(1.7836196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6642163) q[2];
sx q[2];
rz(-1.5216516) q[2];
sx q[2];
rz(1.8159332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3141827) q[1];
sx q[1];
rz(-2.7858509) q[1];
sx q[1];
rz(1.1886503) q[1];
rz(-1.2232699) q[3];
sx q[3];
rz(-1.1545968) q[3];
sx q[3];
rz(-0.89194854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.015335036) q[2];
sx q[2];
rz(-0.33053645) q[2];
sx q[2];
rz(0.34273657) q[2];
rz(0.98894173) q[3];
sx q[3];
rz(-1.9804695) q[3];
sx q[3];
rz(0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2547176) q[0];
sx q[0];
rz(-0.29061341) q[0];
sx q[0];
rz(-1.2683723) q[0];
rz(-2.7791924) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(-1.5836345) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11626205) q[0];
sx q[0];
rz(-2.0206578) q[0];
sx q[0];
rz(2.5309358) q[0];
rz(-1.2940295) q[2];
sx q[2];
rz(-1.74969) q[2];
sx q[2];
rz(-1.8362294) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6889188) q[1];
sx q[1];
rz(-0.59725475) q[1];
sx q[1];
rz(0.026845499) q[1];
rz(-0.5823808) q[3];
sx q[3];
rz(-2.5413168) q[3];
sx q[3];
rz(-1.511661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4871939) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-0.85713345) q[2];
rz(1.0355787) q[3];
sx q[3];
rz(-0.77719694) q[3];
sx q[3];
rz(-2.5100759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421748) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(0.33472043) q[0];
rz(-2.1597629) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(2.2614711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28583254) q[0];
sx q[0];
rz(-1.2783534) q[0];
sx q[0];
rz(0.16181914) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8354957) q[2];
sx q[2];
rz(-1.7405501) q[2];
sx q[2];
rz(1.614326) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0047061027) q[1];
sx q[1];
rz(-2.0396173) q[1];
sx q[1];
rz(-0.15699082) q[1];
rz(-0.91428925) q[3];
sx q[3];
rz(-2.5733272) q[3];
sx q[3];
rz(-1.2340581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52241391) q[2];
sx q[2];
rz(-2.652707) q[2];
sx q[2];
rz(-1.6467113) q[2];
rz(-1.8374247) q[3];
sx q[3];
rz(-0.84916484) q[3];
sx q[3];
rz(-0.18613923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78231597) q[0];
sx q[0];
rz(-0.19141153) q[0];
sx q[0];
rz(3.0707448) q[0];
rz(2.518636) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(-0.43359044) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529652) q[0];
sx q[0];
rz(-0.30634964) q[0];
sx q[0];
rz(3.0393837) q[0];
x q[1];
rz(-1.5135399) q[2];
sx q[2];
rz(-1.861915) q[2];
sx q[2];
rz(-1.8529441) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1939056) q[1];
sx q[1];
rz(-2.8031859) q[1];
sx q[1];
rz(2.848495) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23728063) q[3];
sx q[3];
rz(-10/(7*pi)) q[3];
sx q[3];
rz(-0.90516289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5101461) q[2];
sx q[2];
rz(-2.6084709) q[2];
sx q[2];
rz(0.637429) q[2];
rz(-1.140945) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(-2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7618074) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(-2.8763212) q[0];
rz(-3.1130262) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(2.2144337) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87510159) q[0];
sx q[0];
rz(-1.3463839) q[0];
sx q[0];
rz(1.5988748) q[0];
rz(-pi) q[1];
rz(-1.4248542) q[2];
sx q[2];
rz(-1.2044014) q[2];
sx q[2];
rz(-1.6820358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3520364) q[1];
sx q[1];
rz(-1.4682143) q[1];
sx q[1];
rz(-2.351715) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42334564) q[3];
sx q[3];
rz(-1.4676499) q[3];
sx q[3];
rz(-0.97319095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.022543) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(-1.0724462) q[2];
rz(0.33506814) q[3];
sx q[3];
rz(-3.0266893) q[3];
sx q[3];
rz(-0.2255628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2304614) q[0];
sx q[0];
rz(-1.9879531) q[0];
sx q[0];
rz(-2.352584) q[0];
rz(-2.543653) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(2.423563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58524581) q[0];
sx q[0];
rz(-1.5296773) q[0];
sx q[0];
rz(-1.7264191) q[0];
x q[1];
rz(-1.4903551) q[2];
sx q[2];
rz(-0.1165963) q[2];
sx q[2];
rz(1.2355905) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2788417) q[1];
sx q[1];
rz(-2.2157744) q[1];
sx q[1];
rz(0.74826333) q[1];
rz(0.47649033) q[3];
sx q[3];
rz(-2.9417178) q[3];
sx q[3];
rz(-1.0651922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.614552) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(-1.7238114) q[2];
rz(2.0333911) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(0.38780701) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80477667) q[0];
sx q[0];
rz(-0.62938654) q[0];
sx q[0];
rz(2.8861073) q[0];
rz(-0.77087036) q[1];
sx q[1];
rz(-1.5893385) q[1];
sx q[1];
rz(-1.5482056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324677) q[0];
sx q[0];
rz(-0.45291474) q[0];
sx q[0];
rz(-1.8324773) q[0];
rz(-pi) q[1];
rz(0.86779197) q[2];
sx q[2];
rz(-1.3940587) q[2];
sx q[2];
rz(-0.3294301) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.378506) q[1];
sx q[1];
rz(-2.6623179) q[1];
sx q[1];
rz(2.5162426) q[1];
rz(-pi) q[2];
rz(0.53867619) q[3];
sx q[3];
rz(-2.4823501) q[3];
sx q[3];
rz(-1.074312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-0.74213433) q[2];
sx q[2];
rz(-1.8871657) q[2];
rz(2.6015688) q[3];
sx q[3];
rz(-3.0906782) q[3];
sx q[3];
rz(-2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.222432) q[0];
sx q[0];
rz(-1.5615015) q[0];
sx q[0];
rz(-1.1175565) q[0];
rz(-0.22668214) q[1];
sx q[1];
rz(-0.10348987) q[1];
sx q[1];
rz(-1.6685974) q[1];
rz(-2.1508455) q[2];
sx q[2];
rz(-1.2540091) q[2];
sx q[2];
rz(-2.2768496) q[2];
rz(-1.2516102) q[3];
sx q[3];
rz(-2.2678015) q[3];
sx q[3];
rz(0.24116596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
