OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6731113) q[0];
sx q[0];
rz(-3.0335479) q[0];
sx q[0];
rz(-0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(1.4563814) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1030758) q[0];
sx q[0];
rz(-1.7060192) q[0];
sx q[0];
rz(-1.9413906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26881071) q[2];
sx q[2];
rz(-0.13829002) q[2];
sx q[2];
rz(-2.2594096) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0691245) q[1];
sx q[1];
rz(-2.1241423) q[1];
sx q[1];
rz(-2.8072559) q[1];
rz(-1.637804) q[3];
sx q[3];
rz(-0.57811224) q[3];
sx q[3];
rz(2.1117977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9911389) q[2];
sx q[2];
rz(-3.1296215) q[2];
sx q[2];
rz(-1.0452622) q[2];
rz(2.1446877) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(1.8256942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5528706) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(1.7887822) q[0];
rz(-3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(1.5997684) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38986015) q[0];
sx q[0];
rz(-0.19987389) q[0];
sx q[0];
rz(-0.78011192) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1253042) q[2];
sx q[2];
rz(-1.5920581) q[2];
sx q[2];
rz(2.4580815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5168851) q[1];
sx q[1];
rz(-1.6017645) q[1];
sx q[1];
rz(0.9743486) q[1];
rz(1.0066693) q[3];
sx q[3];
rz(-2.4897277) q[3];
sx q[3];
rz(-1.3019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6359977) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(0.51844281) q[2];
rz(-1.1238267) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(-2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951303) q[0];
sx q[0];
rz(-0.050124425) q[0];
sx q[0];
rz(1.6019524) q[0];
rz(-0.72499544) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(-2.4620788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5121661) q[0];
sx q[0];
rz(-1.1376732) q[0];
sx q[0];
rz(-0.64164759) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1199961) q[2];
sx q[2];
rz(-2.9271002) q[2];
sx q[2];
rz(0.15211764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5522436) q[1];
sx q[1];
rz(-1.2584983) q[1];
sx q[1];
rz(2.8683788) q[1];
x q[2];
rz(1.9539388) q[3];
sx q[3];
rz(-0.62316286) q[3];
sx q[3];
rz(-2.9083657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9069549) q[2];
sx q[2];
rz(-2.8997771) q[2];
sx q[2];
rz(-0.50092906) q[2];
rz(2.6856954) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(-2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86922115) q[0];
sx q[0];
rz(-3.0932194) q[0];
sx q[0];
rz(-0.79917556) q[0];
rz(-0.29798206) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(0.90202773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69660115) q[0];
sx q[0];
rz(-1.7081015) q[0];
sx q[0];
rz(2.3855462) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7219826) q[2];
sx q[2];
rz(-2.107321) q[2];
sx q[2];
rz(-0.46162185) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8405172) q[1];
sx q[1];
rz(-1.4231735) q[1];
sx q[1];
rz(0.0096489659) q[1];
x q[2];
rz(1.6482192) q[3];
sx q[3];
rz(-1.6257878) q[3];
sx q[3];
rz(2.1907326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4681089) q[2];
sx q[2];
rz(-0.031688422) q[2];
sx q[2];
rz(-1.4704977) q[2];
rz(2.9521613) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(-0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162132) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(-2.7498229) q[0];
rz(1.0852934) q[1];
sx q[1];
rz(-3.1317874) q[1];
sx q[1];
rz(-0.36878961) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1802496) q[0];
sx q[0];
rz(-1.215544) q[0];
sx q[0];
rz(0.60549462) q[0];
rz(-pi) q[1];
rz(-2.872345) q[2];
sx q[2];
rz(-1.1675861) q[2];
sx q[2];
rz(2.0079812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49807355) q[1];
sx q[1];
rz(-0.0066702492) q[1];
sx q[1];
rz(1.4590864) q[1];
rz(0.27009876) q[3];
sx q[3];
rz(-0.72805792) q[3];
sx q[3];
rz(-0.61157507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5607249) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(1.3690534) q[2];
rz(0.12618682) q[3];
sx q[3];
rz(-0.43712619) q[3];
sx q[3];
rz(2.0807467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5419902) q[0];
sx q[0];
rz(-0.76102155) q[0];
sx q[0];
rz(1.5916995) q[0];
rz(0.97269336) q[1];
sx q[1];
rz(-2.9284271) q[1];
sx q[1];
rz(2.1125643) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3297076) q[0];
sx q[0];
rz(-2.4548479) q[0];
sx q[0];
rz(2.6558557) q[0];
rz(2.9036361) q[2];
sx q[2];
rz(-1.7081485) q[2];
sx q[2];
rz(2.1417625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29518471) q[1];
sx q[1];
rz(-1.4746397) q[1];
sx q[1];
rz(0.00043649159) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4033086) q[3];
sx q[3];
rz(-1.5351908) q[3];
sx q[3];
rz(-1.6990341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.782393) q[2];
sx q[2];
rz(-1.7155557) q[2];
sx q[2];
rz(0.4314118) q[2];
rz(-2.1443478) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(-0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3417086) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(1.0579911) q[0];
rz(0.82408389) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(0.81738671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1202755) q[0];
sx q[0];
rz(-0.45545721) q[0];
sx q[0];
rz(-1.0767471) q[0];
rz(-3.1233568) q[2];
sx q[2];
rz(-1.5638886) q[2];
sx q[2];
rz(-0.12679174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.83113545) q[1];
sx q[1];
rz(-1.5316446) q[1];
sx q[1];
rz(-1.3292946) q[1];
x q[2];
rz(-2.1002186) q[3];
sx q[3];
rz(-2.2499871) q[3];
sx q[3];
rz(-0.27959945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99007964) q[2];
sx q[2];
rz(-1.9033056) q[2];
sx q[2];
rz(-1.5129169) q[2];
rz(0.40142909) q[3];
sx q[3];
rz(-3.1135058) q[3];
sx q[3];
rz(2.1872971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7734739) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(-0.041944567) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(-1.0036453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40239375) q[0];
sx q[0];
rz(-0.93951997) q[0];
sx q[0];
rz(-1.4556134) q[0];
rz(1.5940158) q[2];
sx q[2];
rz(-2.8041511) q[2];
sx q[2];
rz(2.0793123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81627203) q[1];
sx q[1];
rz(-0.31041103) q[1];
sx q[1];
rz(-2.2903633) q[1];
rz(-2.7549289) q[3];
sx q[3];
rz(-2.4485976) q[3];
sx q[3];
rz(-2.3946144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4369138) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(-1.7436279) q[2];
rz(-2.6513903) q[3];
sx q[3];
rz(-3.0981045) q[3];
sx q[3];
rz(1.1787666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101508) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-0.21139938) q[0];
rz(-2.0231817) q[1];
sx q[1];
rz(-3.1299997) q[1];
sx q[1];
rz(-1.6308019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495962) q[0];
sx q[0];
rz(-1.8798774) q[0];
sx q[0];
rz(-1.2675955) q[0];
rz(-pi) q[1];
rz(-1.1338992) q[2];
sx q[2];
rz(-2.1497576) q[2];
sx q[2];
rz(1.4329662) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76912266) q[1];
sx q[1];
rz(-1.587953) q[1];
sx q[1];
rz(-2.9869798) q[1];
rz(-pi) q[2];
rz(-2.9752067) q[3];
sx q[3];
rz(-1.659044) q[3];
sx q[3];
rz(-1.4180753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3596892) q[2];
sx q[2];
rz(-0.017904559) q[2];
sx q[2];
rz(0.27931279) q[2];
rz(0.8638047) q[3];
sx q[3];
rz(-0.0044718663) q[3];
sx q[3];
rz(2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(1.8260691) q[0];
rz(2.3663991) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(1.388789) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7167146) q[0];
sx q[0];
rz(-1.8535063) q[0];
sx q[0];
rz(1.4155875) q[0];
rz(-pi) q[1];
rz(-2.6493373) q[2];
sx q[2];
rz(-2.6894719) q[2];
sx q[2];
rz(-1.7017572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5701616) q[1];
sx q[1];
rz(-1.5699016) q[1];
sx q[1];
rz(-1.5714297) q[1];
rz(-pi) q[2];
rz(0.41749145) q[3];
sx q[3];
rz(-1.1840828) q[3];
sx q[3];
rz(-1.8698805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3751601) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(0.25894138) q[3];
sx q[3];
rz(-0.18766923) q[3];
sx q[3];
rz(1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5158841) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(-1.6547849) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(1.5776154) q[2];
sx q[2];
rz(-1.7710254) q[2];
sx q[2];
rz(0.22631021) q[2];
rz(0.026080118) q[3];
sx q[3];
rz(-1.9175538) q[3];
sx q[3];
rz(-3.1160311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
