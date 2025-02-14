OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4189897) q[0];
sx q[0];
rz(-0.64282066) q[0];
sx q[0];
rz(1.2195725) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(-2.2396483) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028536) q[0];
sx q[0];
rz(-1.4172638) q[0];
sx q[0];
rz(0.6153591) q[0];
rz(-1.3839533) q[2];
sx q[2];
rz(-1.1964436) q[2];
sx q[2];
rz(0.49546212) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28022596) q[1];
sx q[1];
rz(-1.8609443) q[1];
sx q[1];
rz(-2.9600083) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70848453) q[3];
sx q[3];
rz(-2.2741541) q[3];
sx q[3];
rz(3.0391703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96282643) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(-2.315305) q[2];
rz(1.7360342) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(-2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590782) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(1.4003117) q[0];
rz(2.0179613) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(1.0981015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0535677) q[0];
sx q[0];
rz(-0.19101772) q[0];
sx q[0];
rz(-0.6121401) q[0];
rz(1.2400572) q[2];
sx q[2];
rz(-1.8914701) q[2];
sx q[2];
rz(-1.641524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9713116) q[1];
sx q[1];
rz(-0.14254293) q[1];
sx q[1];
rz(0.55640039) q[1];
rz(-pi) q[2];
rz(-0.16657515) q[3];
sx q[3];
rz(-1.1984571) q[3];
sx q[3];
rz(2.6083715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0565679) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(-1.2497831) q[2];
rz(1.362484) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63610858) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(-2.8662477) q[0];
rz(2.6823726) q[1];
sx q[1];
rz(-2.1870435) q[1];
sx q[1];
rz(-3.088248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5551852) q[0];
sx q[0];
rz(-1.4925028) q[0];
sx q[0];
rz(-0.21348641) q[0];
rz(-pi) q[1];
rz(0.01845999) q[2];
sx q[2];
rz(-1.109913) q[2];
sx q[2];
rz(-0.26319749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1054105) q[1];
sx q[1];
rz(-1.0193258) q[1];
sx q[1];
rz(1.6185733) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1249849) q[3];
sx q[3];
rz(-0.38180581) q[3];
sx q[3];
rz(3.0757381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(1.6374755) q[2];
rz(-1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(0.27568451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212873) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(-2.6570008) q[0];
rz(-1.2424319) q[1];
sx q[1];
rz(-0.74598765) q[1];
sx q[1];
rz(-2.7925083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.718133) q[0];
sx q[0];
rz(-1.0225931) q[0];
sx q[0];
rz(3.0915652) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7789087) q[2];
sx q[2];
rz(-0.35018626) q[2];
sx q[2];
rz(1.9961568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.0047093948) q[1];
sx q[1];
rz(-1.3420958) q[1];
sx q[1];
rz(1.3761702) q[1];
rz(-2.1717668) q[3];
sx q[3];
rz(-1.9221913) q[3];
sx q[3];
rz(1.2343719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7829973) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-0.45004582) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5476371) q[3];
sx q[3];
rz(2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(3.041748) q[0];
rz(-1.3205344) q[1];
sx q[1];
rz(-0.99597275) q[1];
sx q[1];
rz(-2.5340396) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488825) q[0];
sx q[0];
rz(-1.4878055) q[0];
sx q[0];
rz(-1.4957499) q[0];
rz(-pi) q[1];
rz(0.54536087) q[2];
sx q[2];
rz(-1.705745) q[2];
sx q[2];
rz(3.1317186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2098908) q[1];
sx q[1];
rz(-1.7875515) q[1];
sx q[1];
rz(-0.58547975) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0436375) q[3];
sx q[3];
rz(-1.3134911) q[3];
sx q[3];
rz(0.78661455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1943835) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(-0.44140205) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7049103) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(-0.5740903) q[0];
rz(2.2615945) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(-1.7990187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8594151) q[0];
sx q[0];
rz(-2.0382152) q[0];
sx q[0];
rz(0.22478454) q[0];
x q[1];
rz(1.1583352) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(0.5328005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9989661) q[1];
sx q[1];
rz(-1.2104831) q[1];
sx q[1];
rz(0.25288881) q[1];
x q[2];
rz(0.66304147) q[3];
sx q[3];
rz(-1.460307) q[3];
sx q[3];
rz(-0.35943951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(-0.44710844) q[2];
rz(2.1111264) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36987385) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(-2.7287927) q[0];
rz(-1.075047) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(-2.3379751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34016383) q[0];
sx q[0];
rz(-1.0699843) q[0];
sx q[0];
rz(2.1636971) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2718524) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(-0.23245959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10104837) q[1];
sx q[1];
rz(-0.78436995) q[1];
sx q[1];
rz(-0.020009832) q[1];
x q[2];
rz(-1.4945696) q[3];
sx q[3];
rz(-1.7496568) q[3];
sx q[3];
rz(3.0146527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86218086) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(-0.55794445) q[0];
rz(2.5803512) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(-1.7254613) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2673906) q[0];
sx q[0];
rz(-2.6010397) q[0];
sx q[0];
rz(0.99733277) q[0];
rz(-pi) q[1];
rz(-2.2165856) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(-0.031161873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76991612) q[1];
sx q[1];
rz(-1.2425819) q[1];
sx q[1];
rz(-0.10336831) q[1];
x q[2];
rz(-0.21702311) q[3];
sx q[3];
rz(-0.76730928) q[3];
sx q[3];
rz(1.2485152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2531565) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-1.0003132) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(-0.15672556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633453) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-0.42525735) q[1];
sx q[1];
rz(-1.3660627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23200792) q[0];
sx q[0];
rz(-0.9238832) q[0];
sx q[0];
rz(0.034897403) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8810358) q[2];
sx q[2];
rz(-1.2149335) q[2];
sx q[2];
rz(0.17525338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8722265) q[1];
sx q[1];
rz(-0.25372836) q[1];
sx q[1];
rz(2.9969508) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82279295) q[3];
sx q[3];
rz(-1.0446915) q[3];
sx q[3];
rz(1.9578809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1395448) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(0.61919332) q[2];
rz(2.5891417) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904952) q[0];
sx q[0];
rz(-0.18432291) q[0];
sx q[0];
rz(-2.6628394) q[0];
rz(-1.9888196) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(-2.1943888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1789492) q[0];
sx q[0];
rz(-1.5367723) q[0];
sx q[0];
rz(-1.7555139) q[0];
x q[1];
rz(0.2453863) q[2];
sx q[2];
rz(-2.1446501) q[2];
sx q[2];
rz(-1.643744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3403907) q[1];
sx q[1];
rz(-1.2653192) q[1];
sx q[1];
rz(-2.3190772) q[1];
rz(-pi) q[2];
rz(0.60444215) q[3];
sx q[3];
rz(-2.7210208) q[3];
sx q[3];
rz(-0.67422359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08854475) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(-1.0531986) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(-1.5679892) q[2];
sx q[2];
rz(-1.8127828) q[2];
sx q[2];
rz(0.91290963) q[2];
rz(-1.7942747) q[3];
sx q[3];
rz(-2.1383994) q[3];
sx q[3];
rz(2.4169328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
